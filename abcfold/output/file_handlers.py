"""File handlers for different output file types."""

import logging
import re
import warnings
from abc import ABC
from enum import Enum
from pathlib import Path

import gemmi
import numpy as np
import orjson
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue

warnings.filterwarnings("ignore")

logger = logging.getLogger("logger")


class FileTypes(Enum):
    """Enum class for the different file types."""

    NPZ = "npz"
    NPY = "npy"
    CIF = "cif"
    JSON = "json"

    @classmethod
    def values(cls):
        """Get the values of the enum members."""
        return [value.value for value in cls.__members__.values()]


class ModelCount(Enum):
    """Enum class for the different model count types."""

    ALL = "all"
    RESIDUES = "residues"

    @classmethod
    def values(cls):
        """Get the values of the enum members."""
        return [value.value for value in cls.__members__.values()]


class ResidueCountType(Enum):
    """Enum class for the different residue count types."""

    AVERAGE = "average"
    CARBONALPHA = "carbonalpha"
    PHOSPHATE = "phosphate"

    @classmethod
    def values(cls):
        """Get the values of the enum members."""
        return [value.value for value in cls.__members__.values()]


class FileBase(ABC):
    """Abstract base class for the different file types."""

    def __init__(self, pathway: str | Path):
        """File path and suffix names."""
        self.pathway = Path(pathway)
        self.suffix = self.pathway.suffix[1:]

    def __str__(self):
        """Return the string representation of the file path."""
        return str(self.pathway)

    def __repr__(self):
        """Return the class representation of the file path."""
        return f"{self.__class__.__name__}({self.pathway})"


class NpzFile(FileBase):
    """File handler for numpy zip files."""

    def __init__(self, npz_file: str | Path):
        """Object to handle npz files.

        Args:
            npz_file (Union[str, Path]): Path to the npz file

        Attributes:
            npz_file (Path): Path to the npz file
            data (dict): Dictionary containing the data from the npz file

        """
        super().__init__(npz_file)
        self.npz_file = Path(npz_file)
        self.data = self.load_npz_file()

    def load_npz_file(self) -> dict:
        """Load the npz file into a dictionary."""
        return dict(np.load(self.npz_file, allow_pickle=True))


class NpyFile(FileBase):
    """File handler for numpy array files."""

    def __init__(self, npy_file: str | Path):
        """Object to handle npy files.

        Args:
            npy_file (Union[str, Path]): Path to the npy file

        Attributes:
            npy_file (Path): Path to the npy file
            data (np.ndarray): Numpy array containing the data from the np

        """
        super().__init__(npy_file)
        self.npy_file = Path(npy_file)
        self.data = self.load_npy_file()

    def load_npy_file(self) -> np.ndarray:
        """Load the npy file into a numpy array."""
        return np.load(self.npy_file, allow_pickle=True)


class CifFile(FileBase):
    """Cif file handler for structural models."""

    def __init__(self, cif_file: str | Path, input_params: dict | None = None):
        """Object to handle cif files.

        Args:
            cif_file (Union[str, Path]): Path to the cif file
            input_params (Optional[dict]): Dictionary containing the input parameters
            used for the model. This is used to distinguish between ligands and
            sequences

        Attributes:
            cif_file (Path): Path to the cif file
            model (Structure): BioPython structure object containing the model
            atom_plddt_per_chain (dict): Dictionary containing the pLDDT scores for
            each atom
            residue_plddt_per_chain (dict): Dictionary containing the pLDDT scores for
            each residue
            plddts (list): List containing the pLDDT scores for each atom
            residue_plddts (list): List containing the pLDDT scores for each residue
            name (str): Name given to the model

        """
        if input_params is None:
            self.input_params = {}
        else:
            self.input_params = input_params

        super().__init__(cif_file)
        self.cif_file = Path(cif_file)
        self.clashes = 0
        self.clashes_residues = 0
        self.model = self.load_cif_file()
        self.__ligand_plddts = None
        self.__plddts = None
        self.__residue_plddts = None
        self.__h_score = None
        self.__name = self.cif_file.stem

    @property
    def name(self):
        """Filename without suffix."""
        return self.__name

    @name.setter
    def name(self, name: str):
        if not isinstance(name, str):
            logger.error("Name must be a string")
            raise ValueError()
        self.__name = name

    @property
    def plddts(self) -> list[float]:
        """The pLDDT scores for each atom in the model."""
        if self.__plddts is not None:
            return self.__plddts

        self.__plddts = [
            v for plddts in self.get_plddt_per_atom().values() for v in plddts
        ]
        return self.__plddts

    @property
    def residue_plddts(self) -> list[float]:
        """The pLDDT scores for each residue in the model."""
        if self.__residue_plddts is not None:
            return self.__residue_plddts

        self.__residue_plddts = [
            v for plddts in self.get_plddt_per_residue().values() for v in plddts
        ]
        return self.__residue_plddts

    @property
    def average_plddt(self) -> float:
        """The average pLDDT score for the model."""
        return float(np.mean(self.plddts))

    @property
    def ligand_plddts(self) -> dict[str, list[float]]:
        """The pLDDT scores for each ligand in the model."""
        if self.__ligand_plddts is not None:
            return self.__ligand_plddts

        self.__ligand_plddts = self.get_plddt_per_ligand()
        return self.__ligand_plddts

    @property
    def h_score(self) -> int:
        """The H score for the model."""
        if self.__h_score is not None:
            return self.__h_score

        self.__h_score = self.calculate_h_score()
        return self.__h_score

    def load_cif_file(self) -> gemmi.Structure:
        """Load the cif file using Gemmi."""
        st = gemmi.read_structure(str(self.pathway), format=gemmi.CoorFormat.Detect)
        st.setup_entities()
        return st

    @property
    def model1(self) -> gemmi.Model:
        """Get the first model in the structure."""
        return self.model[0]

    def get_chains(self):
        """Alias for self.model1 for compatibility with prev versions of ABCFold2."""
        return self.model1

    def chain_lengths(
        self, mode=ModelCount.RESIDUES, ligand_atoms=False, ptm_atoms=False
    ) -> dict[str, int]:
        """Get the length of each chain in the model.

        Args:
            mode (ModelCount): Enum class specifying the mode to use.
                Note: For ligands the length will always be the number of atoms
            ligand_atoms (bool): Whether to count the number of atoms for ligands.
                If False, the length will be 1 for each ligand chain.
            ptm_atoms (bool): Whether to count the number of atoms for PTMs.
                If False, each PTM will be counted as 1 residue.

        Returns:
            dict: Dictionary containing the chain id and the length of the chain

        Raises:
            ValueError: If the mode is not valid

        """
        model = self.model1
        if mode == ModelCount.ALL or mode == ModelCount.ALL.value:
            chain_lengths: dict = {}
            for chain in model:
                chain: gemmi.Chain
                if chain.name in chain_lengths:
                    chain_lengths[chain.name] += chain.count_atom_sites()
                else:
                    chain_lengths[chain.name] = chain.count_atom_sites()

            return chain_lengths

        elif mode == ModelCount.RESIDUES or mode == ModelCount.RESIDUES.value:
            residue_counts: dict = {}
            for chain in model:
                chain: gemmi.Chain
                chain_polymer = chain.get_polymer()

                # Protein/DNA/RNA chain, count PTMs
                if len(chain_polymer) > 0 and ptm_atoms:
                    residue_counts[chain.name] = sum(
                        1
                        if gemmi.find_tabulated_residue(residue.name).is_standard()
                        else len(residue)
                        for residue in chain
                    )

                elif self.check_ligand(chain):
                    if ligand_atoms:
                        if chain.name in residue_counts:
                            residue_counts[chain.name] += chain.count_atom_sites()
                        else:
                            residue_counts[chain.name] = chain.count_atom_sites()

                    else:
                        residue_counts[chain.name] = 1

                # Protein/DNA/RNA chain, ignore PTMs
                else:
                    residue_counts[chain.name] = len(chain_polymer)

            return residue_counts

        else:
            msg = f"Invalid mode. Please use {', '.join(ModelCount.__members__)}"
            logger.critical(msg)
            raise ValueError()

    def token_residue_ids(self) -> dict:
        """Get the residue ids for each chain in the model.

        Returns:
            dict: Dictionary containing the chain id and the residue ids for each
            chain

        """
        from abcfold.output.utils import flatten

        residue_ids = {}
        for chain in self.model1:
            # In BioPython, residue.id is a tuple of:
            # (hetero flag, sequence identifier, insertion code)
            # In Gemmi, this corresponds to:
            # (res.het_flag, res.seqid.num, res.seqid.icode)
            if self.check_ligand(chain):
                residue_ids[chain.name] = [
                    [residue.seqid.num] for residue in chain for _ in residue
                ]
                continue
            residue_ids[chain.name] = [
                (
                    [residue.seqid.num]
                    if residue.het_flag in {"A", "H"}
                    # if residue.id[0] == " " or residue.id[0] == "H"
                    else [residue.seqid.num for _ in residue]
                )
                for residue in chain
            ]

        residue_ids = {k: flatten(v) for k, v in residue_ids.items()}

        return residue_ids

    def calculate_h_score(self) -> int:
        """Calculate the H score for the model.

        Returns:
            float: The H score for the model

        """
        score = 0
        for i in reversed(range(1, 101)):
            if (100.0 / len(self.plddts)) * np.sum(np.array(self.plddts) >= i) >= i:
                score = i
                break
        return score

    def get_model_sequence_data(self) -> dict[str, str]:
        """Get the sequence for each chain and ligand in the model, used internally for plotting.

        Returns:
            dict : Chain ID and sequence data

        """
        sequence_data = {}
        for chain in self.model1:
            if self.check_ligand(chain):
                # Use the first letter of the atom name as the "sequence" for ligands
                sequence_data[chain.name] = "".join(
                    [atom.name[0] for residue in chain for atom in residue]
                )
            else:
                sequence_data[chain.name] = (
                    chain.get_polymer().make_one_letter_sequence()
                )
        return sequence_data

    def get_plddt_per_atom(self) -> dict[str, list[float]]:
        """Get the pLDDT scores for each atom in the model.

        Returns:
            dict: Dictionary containing the chain id and the pLDDT scores for each atom

        """
        return {
            chain.name: [atom.b_iso for residue in chain for atom in residue]
            for chain in self.model1
        }

    def get_plddt_per_residue(
        self, method=ResidueCountType.AVERAGE.value
    ) -> dict[str, list[float]]:
        """Get the pLDDT scores for each residue in the model.

        Args:
            method (ResidueCountType): Enum class specifying the method to use

        Returns:
            dict: Dictionary containing the chain id and the pLDDT scores for each
            residue

        """
        plddts: dict[str, list[float]] = {}

        if method not in ResidueCountType.values():
            logger.error(
                f"Invalid method. Please use {', '.join(ResidueCountType.__members__)}"
            )
            raise ValueError()

        for chain in self.model1:
            if self.check_ligand(chain):
                # For ligands, return pLDDT per atom
                if chain.name not in plddts:
                    plddts[chain.name] = []
                plddts[chain.name].extend(
                    atom.b_iso for residue in chain for atom in residue
                )

            else:
                for residue in chain:
                    if method == ResidueCountType.AVERAGE.value:
                        score = sum(atom.b_iso for atom in residue) / len(residue)

                    elif method == ResidueCountType.CARBONALPHA.value:
                        score = residue.get_ca().b_iso

                    elif method == ResidueCountType.PHOSPHATE.value:
                        score = residue.get_p().b_iso

                    if chain.name not in plddts:
                        plddts[chain.name] = []
                    plddts[chain.name].append(score)

        plddt_lengths = {k: len(v) for (k, v) in plddts.items()}
        chain_lengths = self.chain_lengths(mode="residues", ligand_atoms=True)
        for chain_id in plddt_lengths:
            assert chain_lengths[chain_id] == plddt_lengths[chain_id], (
                f"{chain_id}, {chain_lengths[chain_id]} != {plddt_lengths[chain_id]}"
            )
        return plddts

    def get_plddt_per_ligand(self) -> dict[str, list[float]]:
        """Get the pLDDT scores for each ligand in the model.

        Returns:
            dict: Dictionary containing the chain id and the pLDDT scores for each atom

        """
        plddt: dict[str, list[float]] = {}
        for chain in self.model1:
            if not self.check_ligand(chain):
                continue

            ligand_plddt = [atom.b_iso for residue in chain for atom in residue]
            if chain.name in plddt:
                plddt[chain.name].extend(ligand_plddt)
            else:
                plddt[chain.name] = ligand_plddt
        return plddt

    def check_ligand(self, chain: gemmi.Chain) -> bool:
        """Check if the chain is a ligand.

        This would have false positives for crystal structures, but for predicted
        structures this should be fine.
        """
        return len(chain.get_ligands()) > 0

    def check_other(self, chain: gemmi.Chain, check_seq_types: list[str]) -> bool:
        """Check if the chain is of a certain type."""
        # TODO: infer sequence types from cif directly
        sequences = self.input_params.get("sequences")
        if sequences is None:
            # logger.warning("Unable to gain sequence information from input file")
            return False
        for sequence in sequences:
            for sequence_type, sequence_data in sequence.items():
                if sequence_type in check_seq_types:
                    if "id" not in sequence_data:
                        continue
                    if hasattr(chain, "name"):
                        chain_id = chain.name
                    else:
                        chain_id = chain
                    if isinstance(sequence_data["id"], str):
                        return chain_id == sequence_data["id"]
                    elif isinstance(sequence_data["id"], list):
                        return chain_id in sequence_data["id"]
        return False

    def relabel_chains(
        self, chain_ids: list[str], link_ids: dict | None = None
    ) -> None:
        """Relabel the chains in the model.

        Args:
            chain_ids (List[str]): List of chain ids to relabel the chains
            link_ids (Optional[dict]): Dictionary containing the chain ids to link
                ligands together. The keys are the new chain ids and the values are
                lists of chain ids to link together.

        Returns:
            None

        """
        # TODO
        chain_ids = chain_ids.copy()
        structure = self.model1
        old_new_chain_id = {}

        if link_ids is None:
            link_ids = {}
        else:
            for new_ids in link_ids.values():
                for new_id in new_ids:
                    chain_ids.pop(chain_ids.index(new_id))

        old_chain_label_counter, new_chain_label_counter = 0, 0

        chain_names = [chain.name for chain in structure]
        while old_chain_label_counter < len(structure):
            chain = chain_ids[new_chain_label_counter]
            old_new_chain_id[chain_names[old_chain_label_counter]] = chain

            # increment the old_chain everytime a chain has been relabelled
            old_chain_label_counter += 1

            if chain in link_ids:
                ligand_no_added = 2
                for _ in link_ids[chain]:
                    old_new_chain_id[chain_names[old_chain_label_counter]] = chain
                    for residue in structure[chain_names[old_chain_label_counter]]:
                        residue.seqid.num = ligand_no_added
                        ligand_no_added += 1
                    old_chain_label_counter += 1

            new_chain_label_counter += 1

        for chain_to_rename in structure:
            chain_to_rename.name = old_new_chain_id[chain_to_rename.name]

        assert old_chain_label_counter == len(self.get_chains()), (
            "Number of chain ids must match the number of chains"
        )
        self.update()

    def update(self):
        """Update the cif file after making changes to the model."""
        self.to_file(self.pathway)
        self = CifFile(self.pathway, self.input_params)

    def reorder_chains(self, new_chain_ids: list[str]):
        """Rearrange the order of chains in a structure model."""
        assert sorted([chain.name for chain in self.model1]) == sorted(new_chain_ids), (
            "The chain ids need to be identical to what is in the model already \
for reordering"
        )
        new_model = gemmi.Model("0")
        for ch in new_chain_ids:
            new_model.add_chain(self.model1[ch])
        new_st = gemmi.Structure()
        new_st.add_model(new_model)
        new_st.meta = self.model.meta
        new_st.raw_remarks = self.model.raw_remarks
        new_st.setup_entities()
        self.model = new_st
        self.update()

    def check_clashes(
        self, threshold: int | float = 3.4, clash_cutoff: float = 0.63
    ) -> tuple[list[tuple[Atom, Atom]], list[tuple[Residue, Residue]]]:
        """Check for clashes between atoms in different chains.

        Args:
            threshold: The distance threshold for a clash.
            clash_cutoff: The cutoff for a clash, as a fraction of the sum of the
                van der Waals radii of the two atoms.

        Returns:
            A list of clashes.

        """
        ns = gemmi.NeighborSearch(
            self.model1, self.model.cell, max_radius=threshold
        ).populate()  # add include_h=False if needed
        cs = gemmi.ContactSearch(threshold)
        cs.ignore = gemmi.ContactSearch.Ignore.SameChain

        # r = a * covalent_radius + b/2, a is multiplier and b is tolerance
        cs.setup_atomic_radii(1.0, 1.5)
        neighbors = cs.find_contacts(ns)

        clashes_atoms, clashes_residues = [], []
        for n in neighbors:
            if n.partner1.atom.name == "C" and n.partner2.atom.name == "N":
                continue
            if n.partner1.atom.name == "N" and n.partner2.atom.name == "C":
                continue
            if (
                n.partner1.atom.name == "SG" and n.partner2.atom.name == "SG"
            ) and n.dist > 1.88:
                continue

            clash_radius = (
                cs.get_radius(n.partner1.atom.element)
                + cs.get_radius(n.partner2.atom.element)
            ) * clash_cutoff
            if n.dist < clash_radius:
                clashes_atoms.append((n.partner1.atom, n.partner2.atom))
                clash_res_pair = (
                    (
                        n.partner1.chain.name,
                        n.partner1.residue.seqid.num,
                        n.partner2.chain.name,
                        n.partner2.residue.seqid.num,
                    )
                    if n.partner1.chain.name < n.partner2.chain.name
                    else (
                        n.partner2.chain.name,
                        n.partner2.residue.seqid.num,
                        n.partner1.chain.name,
                        n.partner1.residue.seqid.num,
                    )
                )
                if clash_res_pair not in clashes_residues:
                    clashes_residues.append(clash_res_pair)

        self.clashes = len(clashes_atoms)
        self.clashes_residues = len(clashes_residues)
        return (clashes_atoms, clashes_residues)

    def get_atoms(self, chain_id=None) -> list[gemmi.Atom]:
        """Get the atoms of the structure."""
        if chain_id is not None:
            return [
                atom
                for chain in self.model1
                for residue in chain
                for atom in residue
                if chain.name == chain_id
            ]
        return [atom for chain in self.model1 for residue in chain for atom in residue]

    def to_file(self, output_file: str | Path) -> None:
        """Save the cif file.

        Args:
            output_file (Union[str, Path]): Path to save the cif file

        Returns:
            None

        """
        out_path = Path(output_file)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        self.model.make_mmcif_document().write_file(str(output_file))
        # io = MMCIFIO()
        # io.set_structure(self.model)

        # # save creates the dictionary
        # io.save(str(output_file))
        # self.__atom_site_label_update(io.dic)
        # self.__ligand_to_hetatm(io.dic)
        # with open(output_file, "w") as f:
        #     io._save_dict(f)

        # self.__single_to_double_quotes(output_file)

    def __single_to_double_quotes(self, file_name: str | Path) -> None:
        new_lines = []
        with open(file_name) as f:
            lines = [line.rstrip() for line in f]

        for line in lines:
            single_quotes = re.compile(r"[']\w+[']{2}")
            new_lines.append(
                single_quotes.sub(lambda x: f'"{x.group()[1:-1]}"', line, count=1)
            )

        with open(file_name, "w") as f:
            f.write("\n".join(new_lines))

    def __atom_site_label_update(self, out_dict):
        atom_site_labels_asym_ids = []
        for chain_id, chain_length in self.chain_lengths(mode="all").items():
            atom_site_labels_asym_ids.extend([chain_id] * chain_length)

        assert len(out_dict["_atom_site.label_asym_id"]) == len(
            atom_site_labels_asym_ids
        ), (
            f"Lengths must be the same, current lengths are \
{len(out_dict['_atom_site.label_asym_id'])} and {len(atom_site_labels_asym_ids)}"
        )

        out_dict["_atom_site.label_asym_id"] = atom_site_labels_asym_ids

        return out_dict

    def __ligand_to_hetatm(self, out_dict):
        atom_site_group_pdb = []
        counter = 0
        for chain_id, chain_length in self.chain_lengths(mode="all").items():
            if self.check_ligand(chain_id):
                atom_site_group_pdb.extend(["HETATM"] * chain_length)
            else:
                atom_site_group_pdb.extend(
                    out_dict["_atom_site.group_PDB"][
                        counter : counter
                        + chain_length  # noqa: E203
                    ]
                )

            counter += chain_length

        assert len(out_dict["_atom_site.group_PDB"]) == len(atom_site_group_pdb)

        out_dict["_atom_site.group_PDB"] = atom_site_group_pdb

        return out_dict


class ConfidenceJsonFile(FileBase):
    """Json file handler for confidence scores."""

    def __init__(self, json_file: str | Path):
        """Object to handle json files.

        Args:
            json_file (Union[str, Path]): Path to the json file

        Attributes:
            json_file (Path): Path to the json file
            data (dict): Dictionary containing the data from the json file

        """
        super().__init__(json_file)
        self.data = self.load_json_file()

    def load_json_file(self):
        """Read the json file into a dict."""
        with open(self.pathway, "rb") as f:
            data = orjson.loads(f.read())

        return data


def superpose_models(
    models_list: list[str | Path], superpose_chains: str | None = None
) -> None:
    """Superpose the models in the list and save them to a new file.

    Args:
        models_list (List[Union[str, Path]]): List of models to superpose
        superpose_chains: Comma-separated chain IDs to use for superposition.
            Defaults to all chains.

    Returns:
        None

    """
    structure = gemmi.read_structure(
        str(models_list[0]), format=gemmi.CoorFormat.Detect
    )
    structure.setup_entities()
    ref_model = structure[0]
    nt_res_names = {"DA", "DT", "DG", "DC", "A", "U", "G", "C", "T"}
    if superpose_chains is not None:
        chain_ids = set(superpose_chains.split(","))
    else:
        chain_ids = {chain.name for chain in ref_model}

    for model in models_list[1:]:
        alt_structure = gemmi.read_structure(str(model), format=gemmi.CoorFormat.Detect)
        alt_structure.setup_entities()
        alt_model = alt_structure[0]

        ref_atoms, alt_atoms = [], []
        for ref_chain, alt_chain in zip(ref_model, alt_model, strict=True):
            if ref_chain.name not in chain_ids:
                continue
            for ref_res, alt_res in zip(ref_chain, alt_chain, strict=True):
                if ref_res.name != alt_res.name or ref_res.seqid != alt_res.seqid:
                    pass

                # Handle nucleotides and proteins differently
                ref_atom_names = {a.name for a in ref_res}
                if ref_res.name in nt_res_names:
                    ref_atoms.append(ref_res.sole_atom("C1'").pos)
                    alt_atoms.append(alt_res.sole_atom("C1'").pos)
                elif "CA" in ref_atom_names:
                    ref_atoms.append(ref_res.get_ca().pos)
                    alt_atoms.append(alt_res.get_ca().pos)
                else:  # Ignore anything else
                    pass

        if len(ref_atoms) == 0 or len(alt_atoms) == 0:
            logger.warning(
                f"No matching atoms found for superposition in {model}. Skipping."
            )
            continue

        sup = gemmi.superpose_positions(ref_atoms, alt_atoms)
        alt_model.transform_pos_and_adp(sup.transform)
        alt_structure.make_mmcif_document().write_file(str(model))
