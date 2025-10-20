"""Schemas for AlphaFold3 input data structures.

https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md
"""

from typing import Annotated

from pydantic import (
    BaseModel,
    NonNegativeInt,
    PlainSerializer,
    PositiveInt,
    model_serializer,
    model_validator,
)


class AF3ProteinModification(BaseModel):
    """Schema for protein modifications."""

    ptmType: str  # CCD code of the PTM
    ptmPosition: PositiveInt  # 1-based index


class AF3NucleotideModification(BaseModel):
    """Schema for DNA/RNA modifications."""

    modificationType: str  # CCD code of the modification
    basePosition: PositiveInt  # 1-based index


class AF3StructuralTemplate(BaseModel):
    """Schema for structural templates.

    Note that the provided mmCIF must contain only a single chain.
    """

    mmcif: str | None = None  # mmCIF string; mutually exclusive with mmcifPath
    mmcifPath: str | None = None
    queryIndices: list[NonNegativeInt]  # 0-based indices in the query sequence
    templateIndices: list[NonNegativeInt]  # 0-based indices in the template sequence

    @model_validator(mode="after")
    def check_mmcif_fields(self):
        """Ensure that exactly one of mmcif or mmcifPath is provided."""
        if (self.mmcif is None) == (self.mmcifPath is None):
            raise ValueError("Exactly one of mmcif or mmcifPath must be provided.")
        return self


class AF3Polymer(BaseModel):
    """Base schema for polymers in AlphaFold3 input."""

    id: str | list[str]  # A, B, ..., Z, AA, BA, CA, ..., ZA, AB, BB, CB, ..., ZB, ...
    sequence: str
    description: str | None = None  # comment describing the chain


class AF3Protein(AF3Polymer):
    """Schema for individual sequences in AlphaFold3 input.

    For recommended MSA options, see:
    https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md#protein-multiple-sequence-alignment
    """

    modifications: list[AF3ProteinModification] | None = None
    # unpairedMsa and unpairedMsaPath are mutually exclusive
    # it's also ok to provide neither. AF3 would error if both are provided
    unpairedMsa: str | None = None  # MSA string in A3M format
    unpairedMsaPath: str | None = (
        None  # Path to MSA file; absolute or relative to the input JSON
    )
    # Same for pairedMsa and pairedMsaPath
    pairedMsa: str | None = None
    pairedMsaPath: str | None = None
    templates: list[AF3StructuralTemplate] | None = None

    @model_validator(mode="after")
    def check_unpaired_msa_fields(self):
        """Ensure that unpairedMsa and unpairedMsaPath are not both provided."""
        if (self.unpairedMsa is not None) and (self.unpairedMsaPath is not None):
            raise ValueError(
                "Only one of unpairedMsa or unpairedMsaPath can be provided."
            )
        return self

    @model_validator(mode="after")
    def check_paired_msa_fields(self):
        """Ensure that pairedMsa and pairedMsaPath are not both provided."""
        if (self.pairedMsa is not None) and (self.pairedMsaPath is not None):
            raise ValueError("Only one of pairedMsa or pairedMsaPath can be provided.")
        return self


class AF3DNA(AF3Polymer):
    """Schema for individual DNA sequences in AlphaFold3 input."""

    modifications: list[AF3NucleotideModification] | None = None


class AF3RNA(AF3Polymer):
    """Schema for individual RNA sequences in AlphaFold3 input.

    For recommended MSA options, see:
    https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md#rna-multiple-sequence-alignment
    """

    modifications: list[AF3NucleotideModification] | None = None
    # unpairedMsa and unpairedMsaPath are mutually exclusive
    # it's also ok to provide neither. AF3 would error if both are provided
    unpairedMsa: str | None = None  # MSA string in A3M format
    unpairedMsaPath: str | None = (
        None  # Path to MSA file; absolute or relative to the input JSON/YAML
    )

    @model_validator(mode="after")
    def check_unpaired_msa_fields(self):
        """Ensure that unpairedMsa and unpairedMsaPath are not both provided."""
        if (self.unpairedMsa is not None) and (self.unpairedMsaPath is not None):
            raise ValueError(
                "Only one of unpairedMsa or unpairedMsaPath can be provided."
            )
        return self


class AF3Ligand(BaseModel):
    """Schema for individual ligand in AlphaFold3 input.

    Ligands can be specified using three formats:

    1. a list of standard CCD codes (e.g., ["HEM", "ZN2"])
    2. a list of custom codes that point to userCCD definitions, e.g. ["LIG-1337"]
    3. a SMILES string that is not in the standard CCD library

    Note that with the SMILES option, you cannot specify covalent bonds to other
    entities as they rely on specific atom names.
    """

    id: str | list[str]  # chain ID(s)
    ccdCodes: list[str] | None = (
        None  # list of standard CCD codes or custom codes pointing to userCCD
    )
    smiles: str | None = None  # optional SMILES string defining the ligand

    description: str | None = None  # comment describing the ligand

    @model_validator(mode="after")
    def check_ccd_smiles_fields(self):
        """Ensure that exactly one of ccdCodes or smiles is provided."""
        if (self.ccdCodes is None) == (self.smiles is None):
            raise ValueError("Exactly one of ccdCodes or smiles must be provided.")
        return self


class AF3Atom(BaseModel):
    """Schema for an atom for specifying bonds."""

    entityId: str  # corresponding to the `id` field for the entity
    residueId: PositiveInt  # 1-based residue index within the chain
    atomName: str  # e.g., "CA", "N", "C", etc.

    @model_serializer
    def serialize_as_list(self) -> list:
        """Serialize as [entityId, resId, atomName]."""
        return [self.entityId, self.residueId, self.atomName]


class AF3BondPair(BaseModel):
    """Schema for bonded atom pairs."""

    atom1: AF3Atom
    atom2: AF3Atom

    @model_serializer
    def serialize_as_list(self) -> list:
        """Serialize as [atom1, atom2]."""
        return [self.atom1, self.atom2]


def ser_af3_sequence(seq: AF3Protein | AF3DNA | AF3RNA) -> dict:
    """Serialize AF3Protein/AF3DNA/AF3RNA as dict with type key."""
    if isinstance(seq, AF3Protein):
        return {"protein": seq}
    elif isinstance(seq, AF3DNA):
        return {"dna": seq}
    elif isinstance(seq, AF3RNA):
        return {"rna": seq}
    else:
        raise TypeError(f"Unsupported sequence type: {type(seq)}")


def ser_af3_sequences(seqs: list[AF3Protein | AF3DNA | AF3RNA]) -> list:
    """Serialize list of AF3Protein/AF3DNA/AF3RNA as list of dicts with type keys."""
    return [ser_af3_sequence(seq) for seq in seqs]


class AF3Input(BaseModel):
    """Input schema for AlphaFold3 structure prediction.

    Refer to: https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md
    """

    name: str
    modelSeeds: list[int]
    sequences: Annotated[
        list[AF3Protein | AF3DNA | AF3RNA], PlainSerializer(ser_af3_sequences)
    ]
    bondedAtomPairs: list[AF3BondPair] | None = None
    # userCCD takes precedence over userCCDPath
    userCCD: str | None = None  # Custom chemical components dictionary
    userCCDPath: str | None = None  # Path to CCD file
    dialect: str = "alphafold3"
    version: PositiveInt = 4
