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

from abcfold.alphafold3.schema import type_key_serializer, type_key_validator
from abcfold.schema import AtomPair as BoltzBond


class BoltzModification(BaseModel):
    """Schema for modifications in Boltz input."""

    position: NonNegativeInt  # 1-based residue index
    ccd: str


class BoltzPolymer(BaseModel):
    """Base schema for polymers in Boltz input.

    DNA/RNA use this schema directly; proteins extend it.
    """

    id: str | list[str]  # A, B, ..., Z, AA, BA, CA, ..., ZA, AB, BB, CB, ..., ZB, ...
    sequence: str
    modifications: list[BoltzModification] | None = None
    cyclic: bool = False


BoltzDNA = BoltzPolymer
BoltzRNA = BoltzPolymer


class BoltzProtein(BoltzPolymer):
    """Schema for protein sequences in Boltz input."""

    # Path to MSA in A3M format; set to "empty" for single-sequence mode
    # Boltz requires an MSA by default
    msa: str | None = None


class BoltzLigand(BaseModel):
    """Schema for ligand sequences in Boltz input."""

    id: str | list[str]  # A, B, ..., Z, AA, BA, CA, ..., ZA, AB, BB, CB, ..., ZB, ...
    smiles: str | None = None
    ccd: str | None = None

    @model_validator(mode="after")
    def check_ccd_smiles_fields(self):
        """Ensure that exactly one of CCD or smiles is provided."""
        if (self.ccd is None) == (self.smiles is None):
            raise ValueError("Exactly one of ccd or smiles must be provided.")
        return self


class BoltzContactResidue(BaseModel):
    """Schema for residues involved in distance constraints in Boltz input."""

    id: str  # chain ID
    residue_idx: NonNegativeInt  # 1-based residue index

    @model_serializer
    def serialize_as_list(self) -> list:
        """Serialize as [chainID, resIdx]."""
        return [self.id, self.residue_idx]


class BoltzContactLigand(BaseModel):
    """Schema for ligand atoms involved in distance constraints in Boltz input."""

    id: str  # chain ID
    atom_name: str

    @model_serializer
    def serialize_as_list(self) -> list:
        """Serialize as [chainID, atomName]."""
        return [self.id, self.atom_name]


def ser_boltz_pocket_contacts(
    contacts: list[BoltzContactResidue | BoltzContactLigand],
) -> list[list[str, int] | list[str, str]]:
    """Serialize list of contact residues/ligands as list of lists."""
    return [contact.model_dump(mode="serialize") for contact in contacts]


class BoltzPocket(BaseModel):
    """Schema for residues associated with binding interaction in Boltz input."""

    binder: str  # chain binding to the pocket
    contacts: Annotated[
        list[BoltzContactResidue | BoltzContactLigand],
        PlainSerializer(ser_boltz_pocket_contacts),
    ]
    max_distance: float = Field(
        default=6.0, ge=4.0, le=20.0
    )  # in Angstroms between any atom in binder and contacts elements
    force: bool = (
        False  # if True, a potential will be used to enforce the pocket constraint
    )


class BoltzContact(BaseModel):
    """Schema for residues associated with binding interaction in Boltz input."""

    token1: BoltzContactResidue | BoltzContactLigand
    token2: BoltzContactResidue | BoltzContactLigand
    max_distance: float = Field(
        default=6.0, ge=4.0, le=20.0
    )  # in Angstroms between any pair of atoms in the two tokens
    force: bool = False


class BoltzStructuralTemplate(BaseModel):
    """Schema for structural templates in Boltz input."""

    cif: str  # path to mmCIF file
    pdb: str | None = None  #
    chain_id: str | list[str] | None = None  # which chain to find a template for
    template_id: str | list[str] | None = None  # which chain in the template to use
    force: bool = False  # if True, a potential will be used to enforce the template
    threshold: float | None = (
        None  # distance (Angstroms) that the prediction can deviate from the template
    )

    @model_validator(mode="after")
    def check_path_fields(self):
        """Ensure that exactly one of cif or pdb is provided."""
        if (self.cif is None) == (self.pdb is None):
            raise ValueError("Exactly one of cif or pdb must be provided.")
        return self


class BoltzAffinity(BaseModel):
    """Schema for affinity prediction in Boltz output."""

    # chain ID of the ligand in the input sequences (with at most 18 atoms)
    # to protein targets. Running against other targets will not raise errors,
    # but results may be unreliable.
    # Running for ligands with >56 atoms is not recommended.
    binder: str


BOLTZ_SEQ_TYPE = {
    "protein": BoltzProtein,
    "dna": BoltzDNA,
    "rna": BoltzRNA,
    "ligand": BoltzLigand,
}


def ser_boltz_sequences(seqs: list[BoltzPolymer | BoltzLigand]) -> list:
    """Serialize list of sequences as list of dicts with type keys."""
    cls2type = {v: k for k, v in BOLTZ_SEQ_TYPE.items()}
    return [type_key_serializer(seq, cls2type) for seq in seqs]


BOLTZ_CONSTRAINT_TYPE = {
    "bond": BoltzBond,
    "pocket": BoltzPocket,
    "contact": BoltzContact,
}


def ser_boltz_constraints(
    constraints: list[BoltzBond | BoltzPocket | BoltzContact] | None,
) -> list | None:
    """Serialize list of constraints as list of dicts with type keys."""
    if constraints is None:
        return None

    cls2type = {v: k for k, v in BOLTZ_CONSTRAINT_TYPE.items()}
    return [type_key_serializer(constraint, cls2type) for constraint in constraints]


BOLTZ_PROPERTY_TYPE = {
    "affinity": BoltzAffinity,
}


def ser_boltz_properties(
    properties: list[BoltzBond | BoltzPocket | BoltzContact] | None,
) -> list | None:
    """Serialize list of properties as list of dicts with type keys."""
    if properties is None:
        return None

    cls2type = {v: k for k, v in BOLTZ_PROPERTY_TYPE.items()}
    return [type_key_serializer(prop, cls2type) for prop in properties]


class BoltzInput(BaseModel):
    """Input schema for Boltz structure prediction.

    Refer to: https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md
    """

    sequences: Annotated[
        list[BoltzProtein | BoltzDNA | BoltzRNA | BoltzLigand],
        PlainSerializer(ser_boltz_sequences),
    ]
    constraints: Annotated[
        list[BoltzBond | BoltzPocket | BoltzContact] | None,
        PlainSerializer(ser_boltz_constraints),
    ] = None
    templates: list[BoltzStructuralTemplate] | None = None
    properties: Annotated[
        list[BoltzAffinity] | None, PlainSerializer(ser_boltz_properties)
    ] = None

    # TODO: add constructor from ABCFold input schema
