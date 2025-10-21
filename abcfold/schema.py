"""Schemas for ABCFold input configs."""

from enum import Enum
from pathlib import Path

from pydantic import (
    BaseModel,
    NonNegativeInt,
    PositiveInt,
    model_serializer,
    model_validator,
)


def type_key_serializer(seq: BaseModel, type_keys: dict[BaseModel, str]):
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
    residue_name: str | None  # Chai requires this for covalent bonds

    @model_serializer
    def serialize_as_list(self) -> list:
        """Serialize as [entityId, resId, atomName]."""
        return [self.chain_id, self.residue_idx, self.atom_name]

    @classmethod
    def init_from_list(cls, atom: list):
        """Initialize from [entityId, resId, atomName]."""
        if len(atom) != 3:
            raise ValueError("Atom list must have exactly 3 elements.")
        return cls(entityId=atom[0], residueId=atom[1], atomName=atom[2])


class AtomPair(BaseModel):
    """Schema for atom pairs."""

    atom1: Atom
    atom2: Atom

    @classmethod
    def init_from_list(cls, atom_pair: list):
        """Initialize from [[entityId, resId, atomName]]."""
        if len(atom_pair) != 2:
            raise ValueError("Atom pair list must have exactly 2 elements.")

        atom1 = Atom.init_from_list(atom_pair[0])
        atom2 = Atom.init_from_list(atom_pair[1])
        return cls(atom1=atom1, atom2=atom2)


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


class ProteinSeq(Polymer):
    """Schema for individual protein sequences."""

    unpaired_msa: str | None = None
    paired_msa: str | None = None
    templates: list[StructuralTemplate] | None = None


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


def load_abcfold_config(conf_file: str) -> ABCFoldConfig:
    """Load ABCFold config from a file."""
    import yaml

    conf_path = Path(conf_file).expanduser().resolve()
    if not conf_path.exists():
        raise FileNotFoundError(f"Config file not found: {conf_path}")
    if conf_path.suffix in {".yml", ".yaml"}:
        with open(conf_path) as f:
            return ABCFoldConfig.model_validate(yaml.safe_load(f))
    elif conf_path.suffix == ".json":
        return ABCFoldConfig.model_validate_json(conf_path.read_bytes())
    else:
        raise ValueError("Unsupported config file format. Use .yaml, .yml, or .json")


def write_config(conf: BaseModel, out_file: str, **kwargs):
    """Write config to a file."""
    import yaml

    yaml.SafeDumper.add_multi_representer(
        Enum,
        yaml.representer.SafeRepresenter.represent_str,
    )

    out_path = Path(out_file).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if "exclude_unset" not in kwargs:
        kwargs["exclude_unset"] = True
    if "exclude_none" not in kwargs:
        kwargs["exclude_none"] = True

    if out_path.suffix in {".yml", ".yaml"}:
        with open(out_path, "w") as f:
            yaml.safe_dump(conf.model_dump(**kwargs), f, sort_keys=False)
    elif out_path.suffix == ".json":
        with open(out_path, "w") as f:
            f.write(conf.model_dump_json(indent=2, **kwargs))
    else:
        raise ValueError("Unsupported config file format. Use .yaml, .yml, or .json")
