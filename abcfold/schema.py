"""Schemas for ABCFold input configs."""

from enum import Enum
from functools import cached_property
from pathlib import Path

from chai_lab.data.dataset.msas.colabfold import generate_colabfold_msas
from chai_lab.data.parsing.msas.aligned_pqt import hash_sequence
from chai_lab.data.parsing.templates.m8 import parse_m8_file
from pydantic import (
    BaseModel,
    NonNegativeInt,
    PositiveInt,
    computed_field,
    model_validator,
)
from tqdm import tqdm

from abcfold.scripts.fetch_rcsb import RCSBDownloader


def type_key_serializer(seq: BaseModel, type_keys: dict[type[BaseModel], str]):
    """Serialize sequences as dict with type key."""
    if type(seq) not in type_keys:
        raise TypeError(f"Unsupported sequence type: {type(seq)}")

    return {type_keys[type(seq)]: seq}


def type_key_validator(key: str, val: BaseModel, type_keys: dict[str, BaseModel]):
    """Validate data based on its type key."""
    if key not in type_keys:
        raise TypeError(f"Unsupported value type: {key}")

    return type_keys[key].model_validate(val)


class Atom(BaseModel):
    """Schema for an atom for specifying bonds."""

    chain_id: str  # corresponding to the `id` field for the entity
    residue_idx: PositiveInt  # 1-based residue index within the chain
    atom_name: str  # e.g., "CA", "N", "C", etc. Follow rdkit for ligands
    residue_name: str | None  # Chai requires this for restraints on proteins

    @classmethod
    def init_from_list(cls, atom: list[str | int]):
        """Initialize from [entityId, resId, atomName]."""
        if len(atom) == 3:
            return cls(chain_id=atom[0], residue_idx=atom[1], atom_name=atom[2])
        elif len(atom) == 4:
            return cls(
                chain_id=atom[0],
                residue_idx=atom[1],
                atom_name=atom[2],
                residue_name=atom[3],
            )
        else:
            raise ValueError(
                "Atom list must have 3 elements and an optional residue_name."
            )


class SequenceModification(BaseModel):
    """Schema for polymer modifications."""

    ccd: str  # CCD code of the PTM
    position: PositiveInt  # 1-based index


class StructuralTemplate(BaseModel):
    """Base schema for structural templates."""

    path: str  # path to the template structure file (mmCIF or PDB)
    # 0-based indices
    query_idx: list[NonNegativeInt] | None = None
    template_idx: list[NonNegativeInt] | None = None
    # IDs in multi-chain templates (not supported in AF3)
    query_chains: list[str] | None = None
    template_chains: list[str] | None = None
    # Boltz-specific fields
    enable_boltz_force: bool = False  # use a potential to enforce the template
    boltz_template_threshold: float | None = (
        None  # distance (Angstroms) that the prediction can deviate from the template
    )


class PolymerType(str, Enum):
    """Enum for polymer types."""

    Protein = "protein"
    DNA = "dna"
    RNA = "rna"

    def __repr__(self):
        """Return string representation when being serialized."""
        return self.value


class Polymer(BaseModel):
    """Base schema for polymers (protein, DNA, and RNA)."""

    seq_type: PolymerType
    id: str | list[str]  # A, B, ..., Z, AA, BA, CA, ..., ZA, AB, BB, CB, ..., ZB, ...
    sequence: str
    modifications: list[SequenceModification] | None = None
    description: str | None = None  # comment describing the chain
    cyclic: bool = False  # Boltz only

    @computed_field
    @cached_property
    def seq_hash(self) -> str:
        """Compute the Chai-style sequence hash."""
        return hash_sequence(self.sequence)


class ProteinSeq(Polymer):
    """Schema for individual protein sequences."""

    msa_dir: str | None = None
    templates: list[StructuralTemplate] | None = None

    @computed_field
    @property
    def unpaired_msa(self) -> str | None:
        """Get path to unpaired MSA file."""
        if self.msa_dir is None:
            return None

        return str(
            Path(self.msa_dir).expanduser().resolve()
            / "a3ms"
            / f"{self.seq_hash}.single.a3m"
        )

    @computed_field
    @property
    def paired_msa(self) -> str | None:
        """Get path to paired MSA file."""
        if self.msa_dir is None:
            return None

        return str(
            Path(self.msa_dir).expanduser().resolve()
            / "a3ms"
            / f"{self.seq_hash}.pair.a3m"
        )


class Ligand(BaseModel):
    """Schema for individual ligands.

    Ligands can be specified using two formats:

    1. a list of standard CCD codes (e.g., ["HEM", "ZN2"])
    2. a SMILES string that is not in the standard CCD library

    Note that with the SMILES option, you cannot specify covalent bonds to other
    entities as they rely on specific atom names.
    """

    id: str | list[str]  # chain ID(s)
    smiles: str | None = None  # optional SMILES string defining the ligand
    ccd: list[str] | None = None  # list of standard CCD codes
    description: str | None = None  # comment describing the ligand

    @model_validator(mode="after")
    def check_ccd_smiles_fields(self):
        """Ensure that exactly one of ccd or smiles is provided."""
        if (self.ccd is None) == (self.smiles is None):
            raise ValueError("Exactly one of ccd or smiles must be provided.")
        return self


class Glycan(BaseModel):
    """Schema for individual glycans."""

    id: str | list[str]  # chain ID(s)
    chai_str: str  # glycan string in Chai notation (modified CCD codes)
    description: str | None = None  # comment describing the glycan


class RestraintType(str, Enum):
    """Enum for restraint types."""

    Covalent = "bond"
    Pocket = "pocket"
    Contact = "contact"

    def __repr__(self):
        """Return string representation when being serialized."""
        return self.value


class Restraint(BaseModel):
    """Schema for distance restraints.

    Note that AF3 only supports bonded restraints.

    In Boltz, the `boltz_binder_chain` should be set to the ligand chain ID that binds
    to the pocket.
    In Chai, the atom that does not belong to `boltz_binder_chain` would be used for
    specifying the pocket, and for the binder chain only the chain ID is needed.
    The atom and residue index information would be ignored.
    Chai also expects the pocket chain to be in Chain B.
    """

    restraint_type: RestraintType
    atom1: Atom
    atom2: Atom
    max_distance: float  # maximum distance (Angstroms); ignored for covalent bonds
    description: str | None = None  # comment describing the restraint

    # Boltz specific fields
    enable_boltz_force: bool = False  # use a potential to enforce the restraint
    boltz_binder_chain: str | None = None  # only used for pocket restraints


class ABCFoldConfig(BaseModel):
    """Config schema for ABCFold."""

    # General settings
    sequences: list[Polymer | ProteinSeq | Ligand | Glycan]
    restraints: list[Restraint] | None = None
    seeds: list[int]

    # Inference parameters
    num_trunk_recycles: int = 3  # Boltz: recycling_steps
    num_diffn_timesteps: int = 200  # Boltz: sampling_steps
    num_diffn_samples: int = 5  # Boltz: diffusion_samples
    num_trunk_samples: int = 1

    # Model-specific settings
    boltz_affinity_binder_chain: str | None = None
    boltz_additional_cli_args: list[str] | None = [
        "--override",
        "--write_full_pae",
        "--write_full_pde",
        # "--use_potentials",
    ]


def load_abcfold_config(conf_file: str | Path) -> ABCFoldConfig:
    """Load ABCFold config from a file."""
    import yaml

    conf_path = Path(conf_file).expanduser().resolve()
    if not conf_path.exists():
        raise FileNotFoundError(f"Config file not found: {conf_path}")
    if conf_path.suffix in {".yml", ".yaml"}:
        with open(conf_path) as f:
            conf = ABCFoldConfig.model_validate(yaml.safe_load(f))
    elif conf_path.suffix == ".json":
        conf = ABCFoldConfig.model_validate_json(conf_path.read_bytes())
    else:
        raise ValueError("Unsupported config file format. Use .yaml, .yml, or .json")

    for i, seq in enumerate(conf.sequences):
        if isinstance(seq, Polymer) and seq.seq_type == PolymerType.Protein:
            conf.sequences[i] = ProteinSeq(**seq.model_dump())

    return conf


def write_config(conf: BaseModel, out_file: str | Path, **kwargs):
    """Write config to a file."""
    import yaml

    yaml.SafeDumper.add_multi_representer(
        Enum,
        yaml.representer.SafeRepresenter.represent_str,
    )

    out_path = Path(out_file).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    for default_arg in ("exclude_unset", "exclude_none", "exclude_computed_fields"):
        if default_arg not in kwargs:
            kwargs[default_arg] = True

    if out_path.suffix in {".yml", ".yaml"}:
        with open(out_path, "w") as f:
            yaml.safe_dump(
                conf.model_dump(**kwargs), f, sort_keys=False, default_flow_style=None
            )
    elif out_path.suffix == ".json":
        with open(out_path, "w") as f:
            f.write(conf.model_dump_json(indent=2, **kwargs))
    else:
        raise ValueError("Unsupported config file format. Use .yaml, .yml, or .json")


def add_msa_to_config(
    conf: ABCFoldConfig,
    out_dir: str | Path,
    chains: set[str] | None = None,
    search_templates: bool = True,
    fetch_templates: bool = True,
    template_cache_dir: Path | None = None,
) -> ABCFoldConfig:
    """Add MSA paths to protein sequences in the config.

    Args:
        conf: ABCFoldConfig object.
        out_dir: Output directory to store MSA files.
        chains: Set of chain IDs to process. If None, process all protein chains.
        search_templates: Whether to search for templates.
        fetch_templates: Whether to fetch mmCIF templates from RCSB.
        template_cache_dir: Directory to cache fetched templates. Defaults to
            ~/.cache/rcsb/ if not set.

    """
    # Figure out which protein sequences to process
    if chains is None:
        chains = {
            c
            for seq in conf.sequences
            if isinstance(seq, ProteinSeq)
            for c in (seq.id if isinstance(seq.id, list) else [seq.id])
        }
    protein_seqs = [
        seq.sequence
        for seq in conf.sequences
        if isinstance(seq, ProteinSeq)
        and any(c in chains for c in (seq.id if isinstance(seq.id, list) else [seq.id]))
    ]

    # Generate MSAs using ColabFold API
    out_path = Path(out_dir).expanduser().resolve()
    out_path.mkdir(parents=True, exist_ok=True)
    _ = generate_colabfold_msas(
        protein_seqs=protein_seqs,
        msa_dir=out_path,
        msa_server_url="https://api.colabfold.com",
        search_templates=search_templates,
        write_a3m_to_msa_dir=True,
    )

    # Fetch templates from RCSB if needed
    if fetch_templates:
        template_cache = (
            Path("~/.cache/rcsb") if template_cache_dir is None else template_cache_dir
        )
        template_cache = template_cache.expanduser().resolve()
        template_cache.mkdir(parents=True, exist_ok=True)
        dl = RCSBDownloader(template_cache, make_subdir=True)

        templates_path = out_path / "all_chain_templates.m8"
        templates_df = parse_m8_file(templates_path)

        template_map: dict[str, list[StructuralTemplate]] = {}
        hash_to_chains: dict[str, list[str]] = {
            seq.seq_hash: (seq.id if isinstance(seq.id, list) else [seq.id])
            for seq in conf.sequences
            if isinstance(seq, ProteinSeq)
        }
        # Enable multithreading here
        for _, r in tqdm(
            templates_df.iterrows(),
            total=templates_df.shape[0],
            desc="Fetching templates",
        ):
            template_pdb_id, template_chain_id = r["subject_id"].split("_")
            template_pdb_id = template_pdb_id.upper()
            template_cif_path = dl.fetch_mmcif(template_pdb_id, file_type="asu")

            if r["query_id"] not in template_map:
                template_map[r["query_id"]] = []
            template_map[r["query_id"]].append(
                StructuralTemplate(
                    path=str(template_cif_path),
                    query_idx=list(range(r["query_start"] - 1, r["query_end"])),
                    template_idx=list(range(r["subject_start"] - 1, r["subject_end"])),
                    query_chains=hash_to_chains[r["query_id"]],
                    template_chains=[template_chain_id],
                )
            )

    # Update config with MSA paths and templates
    (out_path / "templates").mkdir(parents=True, exist_ok=True)
    for i, seq in enumerate(conf.sequences):
        if not isinstance(seq, ProteinSeq):
            continue
        if isinstance(seq.id, str) and seq.id not in chains:
            continue
        elif isinstance(seq.id, list) and not any(c in chains for c in seq.id):
            continue

        conf.sequences[i].msa_dir = str(out_path)

        if seq.seq_hash in template_map:
            custom_templates = seq.templates or []
            conf.sequences[i].templates = custom_templates + template_map[seq.seq_hash]

            # Soft link all templates to out-dir/templates
            for j, t in enumerate(conf.sequences[i].templates):
                template_path = Path(t.path)
                link_path = out_path / "templates" / template_path.name
                if not link_path.exists():
                    link_path.symlink_to(template_path)
                conf.sequences[i].templates[j].path = str(link_path)

    return conf
