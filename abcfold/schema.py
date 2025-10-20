"""Schemas for ABCFold input configs."""

from pydantic import (
    BaseModel,
    PositiveInt,
    model_serializer,
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

    entityId: str  # corresponding to the `id` field for the entity
    residueId: PositiveInt  # 1-based residue index within the chain
    atomName: str  # e.g., "CA", "N", "C", etc.

    @model_serializer
    def serialize_as_list(self) -> list:
        """Serialize as [entityId, resId, atomName]."""
        return [self.entityId, self.residueId, self.atomName]

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
    """Schema for protein modifications."""

    ccd: str  # CCD code of the PTM
    position: PositiveInt  # 1-based index
