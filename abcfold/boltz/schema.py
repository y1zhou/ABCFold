"""Schemas for Boltz2 input YAML.

https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md
"""

from typing import Annotated

from pydantic import (
    BaseModel,
    Field,
    NonNegativeInt,
    PlainSerializer,
    model_serializer,
    model_validator,
)

from abcfold.alphafold3.schema import AF3Atom as Boltz2Atom
from abcfold.alphafold3.schema import type_key_serializer


class Boltz2Modification(BaseModel):
    """Schema for modifications in Boltz2 input."""

    position: NonNegativeInt  # 1-based residue index
    ccd: str


class Boltz2Polymer(BaseModel):
    """Base schema for polymers in Boltz2 input.

    DNA/RNA use this schema directly; proteins extend it.
    """

    id: str | list[str]  # A, B, ..., Z, AA, BA, CA, ..., ZA, AB, BB, CB, ..., ZB, ...
    sequence: str
    modifications: list[Boltz2Modification] | None = None
    cyclic: bool = False


Boltz2DNA = Boltz2Polymer
Boltz2RNA = Boltz2Polymer


class Boltz2Protein(Boltz2Polymer):
    """Schema for protein sequences in Boltz2 input."""

    # Path to MSA in A3M format; set to "empty" for single-sequence mode
    # Boltz2 requires an MSA by default
    msa: str | None = None


class Boltz2Ligand(BaseModel):
    """Schema for ligand sequences in Boltz2 input."""

    id: str | list[str]  # A, B, ..., Z, AA, BA, CA, ..., ZA, AB, BB, CB, ..., ZB, ...
    smiles: str | None = None
    ccd: str | None = None

    @model_validator(mode="after")
    def check_ccd_smiles_fields(self):
        """Ensure that exactly one of CCD or smiles is provided."""
        if (self.ccd is None) == (self.smiles is None):
            raise ValueError("Exactly one of ccd or smiles must be provided.")
        return self


class Boltz2Bond(BaseModel):
    """Schema for covalent bonds in Boltz2 input."""

    atom1: Boltz2Atom
    atom2: Boltz2Atom


class Boltz2ContactResidue(BaseModel):
    """Schema for residues involved in distance constraints in Boltz2 input."""

    id: str  # chain ID
    residue_id: NonNegativeInt  # 1-based residue index

    @model_serializer
    def serialize_as_list(self) -> list:
        """Serialize as [chainID, resIdx]."""
        return [self.id, self.residueId]


class Boltz2ContactLigand(BaseModel):
    """Schema for ligand atoms involved in distance constraints in Boltz2 input."""

    id: str  # chain ID
    atom_name: str

    @model_serializer
    def serialize_as_list(self) -> list:
        """Serialize as [chainID, atomName]."""
        return [self.id, self.atom_name]


def ser_boltz2_pocket_contacts(
    contacts: list[Boltz2ContactResidue | Boltz2ContactLigand],
) -> list[list[str, int] | list[str, str]]:
    """Serialize list of contact residues/ligands as list of lists."""
    return [contact.model_dump(mode="serialize") for contact in contacts]


class Boltz2Pocket(BaseModel):
    """Schema for residues associated with binding interaction in Boltz2 input."""

    binder: str  # chain binding to the pocket
    contacts: Annotated[
        list[Boltz2ContactResidue | Boltz2ContactLigand],
        PlainSerializer(ser_boltz2_pocket_contacts),
    ]
    max_distance: float = Field(
        default=6.0, ge=4.0, le=20.0
    )  # in Angstroms between any atom in binder and contacts elements
    force: bool = (
        False  # if True, a potential will be used to enforce the pocket constraint
    )


class Boltz2Contact(BaseModel):
    """Schema for residues associated with binding interaction in Boltz2 input."""

    token1: Boltz2ContactResidue | Boltz2ContactLigand
    token2: Boltz2ContactResidue | Boltz2ContactLigand
    max_distance: float = Field(
        default=6.0, ge=4.0, le=20.0
    )  # in Angstroms between any pair of atoms in the two tokens
    force: bool = False


class Boltz2StructuralTemplate(BaseModel):
    """Schema for structural templates in Boltz2 input."""

    cif: str  # path to mmCIF file
    pdb: str | None = None  #
    chain_id: str | list[str] | None = None  # which chain to find a template for
    template_id: str | list[str] | None = None  # which chain in the template to use
    force: bool = False  # if True, a potential will be used to enforce the template
    threshold: float | None = (
        None  # distance (in Angstroms) that the prediction can deviate from the template
    )

    @model_validator(mode="after")
    def check_path_fields(self):
        """Ensure that exactly one of cif or pdb is provided."""
        if (self.cif is None) == (self.pdb is None):
            raise ValueError("Exactly one of cif or pdb must be provided.")
        return self


class Boltz2Affinity(BaseModel):
    """Schema for affinity prediction in Boltz2 output."""

    # chain ID of the ligand in the input sequences (with at most 128 atoms)
    # to protein targets. Running against other targets will not raise errors,
    # but results may be unreliable.
    # Running for ligands with >56 atoms is not recommended.
    binder: str


def ser_boltz2_sequences(seqs: list[Boltz2Polymer | Boltz2Ligand]) -> list:
    """Serialize list of sequences as list of dicts with type keys."""
    type_keys = {
        Boltz2Protein: "protein",
        Boltz2DNA: "dna",
        Boltz2RNA: "rna",
        Boltz2Ligand: "ligand",
    }
    return [type_key_serializer(seq, type_keys) for seq in seqs]


def ser_boltz2_constraints(
    constraints: list[Boltz2Bond | Boltz2Pocket | Boltz2Contact] | None,
) -> list | None:
    """Serialize list of constraints as list of dicts with type keys."""
    if constraints is None:
        return None

    type_keys = {
        Boltz2Bond: "bond",
        Boltz2Pocket: "pocket",
        Boltz2Contact: "contact",
    }
    return [type_key_serializer(constraint, type_keys) for constraint in constraints]


def ser_boltz2_properties(
    properties: list[Boltz2Bond | Boltz2Pocket | Boltz2Contact] | None,
) -> list | None:
    """Serialize list of properties as list of dicts with type keys."""
    if properties is None:
        return None

    type_keys = {
        Boltz2Affinity: "affinity",
    }
    return [type_key_serializer(prop, type_keys) for prop in properties]


class Boltz2Input(BaseModel):
    """Input schema for Boltz2 structure prediction.

    Refer to: https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md
    """

    sequences: Annotated[
        list[Boltz2Protein | Boltz2DNA | Boltz2RNA | Boltz2Ligand],
        PlainSerializer(ser_boltz2_sequences),
    ]
    constraints: Annotated[
        list[Boltz2Bond | Boltz2Pocket | Boltz2Contact] | None,
        PlainSerializer(ser_boltz2_sequences),
    ] = None
    templates: list[Boltz2StructuralTemplate] | None = None
    properties: Annotated[
        list[Boltz2Affinity] | None, PlainSerializer(ser_boltz2_properties)
    ] = None

    # TODO: add constructor from ABCFold input schema
