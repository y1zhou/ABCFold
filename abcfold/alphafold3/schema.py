"""Schemas for AlphaFold3 input data structures.

https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md
"""

from pydantic import BaseModel


class AF3ProteinModification(BaseModel):
    """Schema for protein modifications."""

    ptmType: str  # CCD code of the PTM
    ptmPosition: int  # 1-based index


class AF3NucleotideModification(BaseModel):
    """Schema for DNA/RNA modifications."""

    modificationType: str  # CCD code of the modification
    basePosition: int  # 1-based index


class AF3StructuralTemplate(BaseModel):
    """Schema for structural templates.

    Note that the provided mmCIF must contain only a single chain.
    """

    # TODO: check for mutually exclusive mmcif and mmcifPath
    mmcif: str | None = None  # mmCIF string; mutually exclusive with mmcifPath
    mmcifPath: str | None = None
    queryIndices: list[int]  # 0-based indices in the query sequence
    templateIndices: list[int]  # 0-based indices in the template sequence


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
    # it's also ok to provide neither
    # AF3 would error if both are provided; here unpairedMsa takes precedence
    unpairedMsa: str | None = None  # MSA string in A3M format
    unpairedMsaPath: str | None = (
        None  # Path to MSA file; absolute or relative to the input JSON/YAML
    )
    # Same for pairedMsa and pairedMsaPath
    pairedMsa: str | None = None
    pairedMsaPath: str | None = None
    templates: list[AF3StructuralTemplate] | None = None


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
    # it's also ok to provide neither
    # AF3 would error if both are provided; here unpairedMsa takes precedence
    unpairedMsa: str | None = None  # MSA string in A3M format
    unpairedMsaPath: str | None = (
        None  # Path to MSA file; absolute or relative to the input JSON/YAML
    )


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
    # TODO: check for mutually exclusive ccdCodes and smiles
    ccdCodes: list[str] | None = (
        None  # list of standard CCD codes or custom codes pointing to userCCD
    )
    smiles: str | None = None  # optional SMILES string defining the ligand

    description: str | None = None  # comment describing the ligand


class AF3Atom(BaseModel):
    """Schema for an atom for specifying bonds."""

    entityId: str  # corresponding to the `id` field for the entity
    residueId: int  # 1-based residue index within the chain
    atomName: str  # e.g., "CA", "N", "C", etc.


class AF3BondPair(BaseModel):
    """Schema for bonded atom pairs."""

    # TODO: serialize as [[entityId, resId, atomName], [entityId, resId, atomName]]
    atom1: AF3Atom
    atom2: AF3Atom


class AF3Input(BaseModel):
    """Input schema for AlphaFold3 structure prediction.

    Refer to: https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md
    """

    name: str
    modelSeeds: list[int]
    # TODO: serialize as list of {"protein|dna|rna": AF3Protein|AF3DNA|AF3RNA}
    sequences: list[AF3Protein | AF3DNA | AF3RNA]
    bondedAtomPairs: list[AF3BondPair] | None = None
    # userCCD takes precedence over userCCDPath
    userCCD: str | None = None  # Custom chemical components dictionary
    userCCDPath: str | None = None  # Path to CCD file
    dialect: str = "alphafold3"
    version: int = 4
