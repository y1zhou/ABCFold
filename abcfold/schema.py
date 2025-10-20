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
    atom_name: str  # e.g., "CA", "N", "C", etc.

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


class RestraintType(str, Enum):
    """Enum for restraint types."""

    Covalent = "bond"
    Pocket = "pocket"
    Contact = "contact"


class RestraintAtom(BaseModel):
    """Schema for an atom in a restraint."""

    chain_id: str  # corresponding to the `id` field for the entity
    residue_name: str | None
    residue_idx: PositiveInt | None  # 1-based residue index within the chain
    atom_name: str  # e.g., "CA", "N", "C", etc. Follow rdkit for ligands


class Restraint(BaseModel):
    """Schema for distance restraints.

    Note that AF3 only supports bonded restraints.
    """

    restraint_type: RestraintType
    atom1: RestraintAtom
    atom2: RestraintAtom
    max_distance: float  # maximum distance (Angstroms); ignored for covalent bonds
    description: str | None = None  # comment describing the restraint

    # Boltz specific fields
    enable_boltz_force: bool = False  # use a potential to enforce the restraint
    boltz_binder_chain: str | None = None


class ABCFoldConfig(BaseModel):
    """Config schema for ABCFold."""

    # General settings
    sequences: list[Polymer | Ligand]
    bonds: list[AtomPair] | None = None
    restraints: list[Restraint] | None = None
    seeds: list[int]

    # Boltz-specific settings
    boltz_affinity_binder_chain: str | None = None
    additional_boltz_cli_args: list[str] | None = None


def load_abcfold_config(conf_file: str) -> ABCFoldConfig:
    """Load ABCFold config from a file."""
    import json

    import yaml

    conf_path = Path(conf_file).expanduser().resolve()
    if not conf_path.exists():
        raise FileNotFoundError(f"Config file not found: {conf_path}")
    if conf_path.suffix in {".yml", ".yaml"}:
        with open(conf_path) as f:
            conf_dict = yaml.safe_load(f)
    elif conf_path.suffix == ".json":
        with open(conf_path) as f:
            conf_dict = json.load(f)
    else:
        raise ValueError("Unsupported config file format. Use .yaml, .yml, or .json")

    return ABCFoldConfig.model_validate(conf_dict)


def write_config(conf: BaseModel, out_file: str, **kwargs):
    """Write config to a file."""
    import json

    import yaml

    out_path = Path(out_file).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    conf_dict = conf.model_dump(**kwargs)
    if out_path.suffix in {".yml", ".yaml"}:
        with open(out_path, "w") as f:
            yaml.safe_dump(conf_dict, f)
    elif out_path.suffix == ".json":
        with open(out_path, "w") as f:
            json.dump(conf_dict, f, indent=2)
    else:
        raise ValueError("Unsupported config file format. Use .yaml, .yml, or .json")
