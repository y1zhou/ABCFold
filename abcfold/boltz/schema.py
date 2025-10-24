"""Schemas for Boltz2 input YAML.

https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md
"""

from pathlib import Path
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
from abcfold.schema import (
    ABCFoldConfig,
    Atom,
    Glycan,
    Ligand,
    Polymer,
    PolymerType,
    ProteinSeq,
    RestraintType,
)
from abcfold.schema import Ligand as BoltzLigand
from abcfold.schema import SequenceModification as BoltzModification


class BoltzAtom(Atom):
    """Schema for atoms in Boltz input.

    Extends Atom with serializer.
    """

    @model_serializer
    def serialize_as_list(self) -> list:
        """Serialize as [entityId, resId, atomName]."""
        return [self.chain_id, self.residue_idx, self.atom_name]


class BoltzBond(BaseModel):
    """Schema for atom pairs."""

    atom1: BoltzAtom
    atom2: BoltzAtom

    @classmethod
    def init_from_list(cls, atom_pair: list[list[str | int]]):
        """Initialize from [[entityId, resId, atomName]]."""
        if len(atom_pair) != 2:
            raise ValueError("Atom pair list must have exactly 2 elements.")

        atom1 = BoltzAtom.init_from_list(atom_pair[0])
        atom2 = BoltzAtom.init_from_list(atom_pair[1])
        return cls(atom1=atom1, atom2=atom2)

    @classmethod
    def init_from_atoms(cls, atom1: Atom, atom2: Atom):
        """Initialize from [[entityId, resId, atomName]]."""
        dump_fields = ("chain_id", "residue_idx", "atom_name")
        atom1 = BoltzAtom(**atom1.model_dump(include=dump_fields))
        atom2 = BoltzAtom(**atom2.model_dump(include=dump_fields))
        return cls(atom1=atom1, atom2=atom2)


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
) -> list[list[str | int] | list[str]]:
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

    cif: str | None = None  # path to mmCIF file
    pdb: str | None = None
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

    @model_validator(mode="after")
    def check_force_fields(self):
        """Ensure that threshold is set when force==True."""
        if self.force and self.threshold is None:
            raise ValueError("threshold must be set when force is set to True.")

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


def ser_boltz_properties(properties: list[BoltzAffinity] | None) -> list | None:
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

    @classmethod
    def init_from_boltz_yaml(cls, dat: dict):
        """Initialize from Boltz input YAML dictionary."""
        seqs = []
        seq_types = {}
        for seq in dat["sequences"]:
            for k, v in seq.items():
                seq_model = type_key_validator(k, v, BOLTZ_SEQ_TYPE)
                seqs.append(seq_model)
                if isinstance(seq_model.id, list):
                    for chain_id in seq_model.id:
                        seq_types[chain_id] = k
                else:
                    seq_types[seq_model.id] = k

        dat["sequences"] = seqs

        constraints = []

        def _build_contact_elem_from_list(chain: str, elem_id: str | int):
            if seq_types[chain] == "ligand":
                return BoltzContactLigand(id=chain, atom_name=elem_id)
            else:
                return BoltzContactResidue(id=chain, residue_idx=elem_id)

        for constraint in dat.get("constraints", []):
            for k, v in constraint.items():
                if k == "bond":
                    constraints.append(BoltzBond.init_from_list(v))
                elif k == "pocket":
                    contact_elems = [
                        _build_contact_elem_from_list(*elem) for elem in v["contacts"]
                    ]
                    v["contacts"] = contact_elems
                    constraints.append(BoltzPocket.model_validate(v))
                elif k == "contact":
                    v["token1"] = _build_contact_elem_from_list(*v["token1"])
                    v["token2"] = _build_contact_elem_from_list(*v["token2"])
                    constraints.append(BoltzContact.model_validate(v))
                else:
                    raise ValueError(f"Unknown constraint type: {k}")
        if constraints:
            dat["constraints"] = constraints

        props = []
        for prop in dat.get("properties", []):
            for k, v in prop.items():
                props.append(type_key_validator(k, v, BOLTZ_PROPERTY_TYPE))
        if props:
            dat["properties"] = props

        return cls.model_validate(dat)


def merge_chai_msa_to_csv(
    unpaired_msa_file: str | Path | None,
    paired_msa_file: str | Path | None,
    msa_id: str,
    out_dir: str | Path,
) -> Path:
    """Merge unpaired and paired MSAs into a single CSV file for Boltz.

    Adapted from Boltz's own MSA processing code: boltz.main.compute_msa

    Args:
        unpaired_msa_file: Path to unpaired MSA file in A3M format.
        paired_msa_file: Path to paired MSA file in A3M format.
        msa_id: Output file name. Boltz uses `{target_id}_{entity_id}.csv`.
        out_dir: Directory to save the output CSV file.

    """
    out_dir_path = Path(out_dir).expanduser().resolve()
    out_dir_path.mkdir(parents=True, exist_ok=True)
    out_file = out_dir_path / f"{msa_id}.csv"

    # fast fail if unpaired MSA file does not exist
    if unpaired_msa_file is None:
        raise ValueError("Unpaired MSA file must be provided.")
    unpaired_path = Path(unpaired_msa_file).expanduser().resolve()
    if not unpaired_path.exists():
        raise FileNotFoundError(f"Unpaired MSA file not found: {unpaired_path}")

    if paired_msa_file is not None:
        paired_path = Path(paired_msa_file).expanduser().resolve()
        if not paired_path.exists():
            raise FileNotFoundError(f"Paired MSA file not found: {paired_path}")
        with open(paired_path) as f:
            paired = f.read().strip().splitlines()

        # ignore headers
        # Boltz also does subsampling here but we skip that for now
        paired = paired[1::2]

        # set key per row and remove empty sequences
        keys = [idx for idx, s in enumerate(paired) if s != "-" * len(s)]
        paired = [s for s in paired if s != "-" * len(s)]
    else:
        paired = []
        keys = []

    # combine paired-unpaired sequences
    with open(unpaired_path) as f:
        unpaired = f.read().strip().splitlines()

    unpaired = unpaired[1::2]
    if paired:
        # ignore query seq
        unpaired = unpaired[1:]

    combined_seqs = paired + unpaired
    combined_keys = keys + [-1] * len(unpaired)

    # Dump to CSV
    csv_str = ["key,sequence"] + [
        f"{key},{seq}" for key, seq in zip(combined_keys, combined_seqs, strict=True)
    ]
    with open(out_file, "w") as f:
        f.write("\n".join(csv_str))

    return out_file


def abcfold_to_boltz(conf: ABCFoldConfig, msa_dir: str | Path) -> BoltzInput:
    """Convert ABCFold config to Boltz input schema.

    Args:
        conf: ABCFold configuration.
        msa_dir: Directory to save MSA CSV files for Boltz.

    """
    # Sequences
    seqs: list[BoltzProtein | BoltzDNA | BoltzRNA | BoltzLigand] = []
    seq_types: dict[str, str] = {}

    # Boltz stores templates in a separate field
    # TODO: Chai loads at most 4 templates per chain; Boltz has no such limit but
    # excessive templates may overwhelm the GPU memory
    templates: list[BoltzStructuralTemplate] = []

    msa_dir_path = Path(msa_dir).expanduser().resolve()
    msa_dir_path.mkdir(parents=True, exist_ok=True)
    for seq in conf.sequences:
        # Record sequence types for each input chain
        if isinstance(seq, Glycan):
            raise ValueError("Glycans are not yet supported in Boltz input.")

        if isinstance(seq.id, list):
            for chain_id in seq.id:
                seq_types[chain_id] = (
                    seq.seq_type.value if isinstance(seq, Polymer) else "ligand"
                )
        else:
            seq_types[seq.id] = (
                seq.seq_type.value if isinstance(seq, Polymer) else "ligand"
            )

        # Map ABCFold sequences to Boltz sequences
        if isinstance(seq, Ligand):
            seqs.append(seq)
        elif isinstance(seq, ProteinSeq):
            if seq.msa_dir is not None:
                msa_csv_path = merge_chai_msa_to_csv(
                    seq.unpaired_msa,
                    seq.paired_msa,
                    msa_id=seq.seq_hash,
                    out_dir=msa_dir_path,
                )
            else:
                msa_csv_path = "empty"
            seqs.append(
                BoltzProtein(
                    id=seq.id,
                    sequence=seq.sequence,
                    modifications=seq.modifications,
                    cyclic=seq.cyclic,
                    msa=str(msa_csv_path),
                )
            )

            # Templates
            if not seq.templates:
                continue
            for tmpl in seq.templates:
                tmpl_path = Path(tmpl.path).expanduser().resolve()
                cif_path, pdb_path = None, None
                match "".join(tmpl_path.suffixes[-2:]).lower():
                    case ".cif" | ".cif.gz":
                        cif_path = str(tmpl_path)
                    case ".pdb" | ".pdb.gz":
                        pdb_path = str(tmpl_path)
                    case _:
                        raise ValueError(
                            f"Unsupported template file format: {tmpl_path}"
                        )
                templates.append(
                    # TODO: Boltz uses gemmi.structure.entities instead of subchains
                    # See boltz.data.parse.mmcif.parse_mmcif
                    BoltzStructuralTemplate(
                        cif=cif_path,
                        pdb=pdb_path,
                        chain_id=tmpl.query_chains,
                        # template_id=tmpl.template_chains,
                        force=tmpl.enable_boltz_force,
                        threshold=tmpl.boltz_template_threshold,
                    )
                )
        elif seq.seq_type == PolymerType.DNA:
            seqs.append(
                BoltzDNA(
                    **seq.model_dump(
                        include={"id", "sequence", "modifications", "cyclic"}
                    )
                )
            )
        elif seq.seq_type == PolymerType.RNA:
            seqs.append(
                BoltzRNA(
                    **seq.model_dump(
                        include={"id", "sequence", "modifications", "cyclic"}
                    )
                )
            )
        else:
            # Glycans not yet supported in Boltz
            raise ValueError(f"Unsupported polymer type: {seq.seq_type}")

    # Constraints
    def _build_contact_elem_from_list(elem: Atom):
        if seq_types[elem.chain_id] == "ligand":
            return BoltzContactLigand(id=elem.chain_id, atom_name=elem.atom_name)
        else:
            return BoltzContactResidue(id=elem.chain_id, residue_idx=elem.residue_idx)

    constraints = []
    for restraint in conf.restraints or []:
        if restraint.restraint_type == RestraintType.Covalent:
            constraints.append(
                BoltzBond.init_from_atoms(atom1=restraint.atom1, atom2=restraint.atom2)
            )
        elif restraint.restraint_type == RestraintType.Pocket:
            if restraint.boltz_binder_chain is None:
                raise ValueError(
                    f"boltz_binder_chain must be specified for pocket restraints: {restraint}"
                )
            constraints.append(
                BoltzPocket(
                    binder=restraint.boltz_binder_chain,
                    contacts=[
                        _build_contact_elem_from_list(restraint.atom1),
                        _build_contact_elem_from_list(restraint.atom2),
                    ],
                    max_distance=restraint.max_distance,
                    force=restraint.enable_boltz_force,
                )
            )
        elif restraint.restraint_type == RestraintType.Contact:
            constraints.append(
                BoltzContact(
                    token1=_build_contact_elem_from_list(restraint.atom1),
                    token2=_build_contact_elem_from_list(restraint.atom2),
                    max_distance=restraint.max_distance,
                    force=restraint.enable_boltz_force,
                )
            )
        else:
            raise ValueError(f"Unsupported restraint type: {restraint.restraint_type}")

    # Properties
    properties = []
    if conf.boltz_affinity_binder_chain is not None:
        properties.append(BoltzAffinity(binder=conf.boltz_affinity_binder_chain))

    boltz_input = BoltzInput(
        sequences=seqs,
        constraints=constraints if constraints else None,
        templates=templates if templates else None,
        properties=properties if properties else None,
    )

    return boltz_input
