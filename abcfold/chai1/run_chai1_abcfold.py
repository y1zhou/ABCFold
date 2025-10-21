"""Run Chai-1 using ABCFold configuration file."""

import logging
from pathlib import Path

import pandas as pd

from abcfold.schema import (
    ABCFoldConfig,
    Atom,
    Glycan,
    Ligand,
    Polymer,
    RestraintType,
    SequenceModification,
    load_abcfold_config,
)

logger = logging.getLogger("logger")


class ChaiConfig:
    """Builder for Chai-1 config files from ABCFold config."""

    def __init__(
        self,
        conf: ABCFoldConfig,
        workdir: str | Path,
        run_id: str,
        ccd_lib_dir: str | Path | None = None,
    ):
        """Initialize ChaiConfig builder.

        Args:
            conf (ABCFoldConfig): ABCFold configuration object.
            workdir (str | Path): Output directory for Chai-1 inputs.
            run_id (str): Unique identifier for the run.
            ccd_lib_dir: Boltz cache directory for CCD ligands.

        """
        self.conf = conf
        self.workdir = Path(workdir).expanduser().resolve()
        self.run_id = run_id
        self.ccd_lib_dir = (
            Path(ccd_lib_dir).expanduser().resolve()
            if ccd_lib_dir is not None
            else None
        )

        self.workdir.mkdir(parents=True, exist_ok=True)

    def generate_chai_inputs(self) -> dict[str, tuple[str, str]]:
        """Generate all Chai-1 input files."""
        # File paths
        self.fasta = self.workdir / f"{self.run_id}.fasta"
        self.restraints = self.workdir / f"{self.run_id}.restraints"

        self.seq_metadata: dict[str, tuple[str, str]] = self.generate_chai_fasta(
            self.fasta, self.ccd_lib_dir
        )
        self.generate_chai_restraints()

    def generate_chai_fasta(
        self, ccd_lib_dir: str | Path | None = None
    ) -> dict[str, tuple[str, str]]:
        """Build Chai-1 style FASTA from ABCFold config.

        `ccd_lib_dir` is needed if their are CCD ligands in the input.
        Chai only supports SMILES ligands.

        Note that Chai-1 outputs re-assigns chain IDs based on the order of sequences
        in the FASTA file. So the chain IDs in the ABCFold config may not be preserved.

        Returns:
            Mapping from original chain IDs to Chai-assigned chain IDs and entity types:
                {original_chain_id: (chai_chain_id, entity_type)}

        """
        entries: list[str] = []
        for seq in self.conf.sequences:
            if isinstance(seq, Glycan):
                entries.extend(self._build_fasta_entry("glycan", seq.id, seq.chai_str))

            elif isinstance(seq, Ligand):
                if seq.ccd:
                    if ccd_lib_dir is None:
                        raise ValueError(
                            "ccd_lib_dir must be provided when CCD ligands are used."
                        )
                    for chain_id, ccd_code in zip(seq.id, seq.ccd, strict=True):
                        smiles = self.ccd_to_smiles(ccd_code, ccd_lib_dir)
                        entries.extend(
                            self._build_fasta_entry("ligand", chain_id, smiles)
                        )
                else:
                    entries.extend(
                        self._build_fasta_entry("ligand", seq.id, seq.smiles)
                    )

            elif isinstance(seq, Polymer):
                entries.extend(
                    self._build_fasta_entry(
                        seq.seq_type.value,
                        seq.id,
                        self._add_modifications_to_sequence(
                            seq.sequence, seq.modifications
                        ),
                    )
                )
            else:
                raise ValueError(f"Unsupported sequence type: {type(seq)}")

        with open(self.fasta, "w") as f:
            f.write("\n".join(entries) + "\n")

        chai_chain_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        seq_metadata = {}
        for i, entry in enumerate(entries):
            entity_type, original_chain_id = entry.splitlines()[0].split("|")
            seq_metadata[original_chain_id] = (
                chai_chain_ids[i],
                entity_type[1:],  # remove '>' from seq_type
            )

        return seq_metadata

    def generate_chai_restraints(self):
        """Generate Chai-1 style restraints CSV from ABCFold config."""
        if not self.conf.restraints:
            self.restraints = None
            return

        restraint_cols = [
            "restraint_id",
            "chainA",
            "res_idxA",
            "chainB",
            "res_idxB",
            "connection_type",
            "confidence",  # ignored
            "min_distance_angstrom",  # ignored
            "max_distance_angstrom",
            "comment",  # ignored
        ]
        entries = []
        for r in self.conf.restraints:
            chain1_id, chain1_type = self.seq_metadata[r.atom1.chain_id]
            chain2_id, chain2_type = self.seq_metadata[r.atom2.chain_id]

            if chain1_id == chain2_id:
                raise ValueError(f"Restraints must be inter-chain for Chai-1, got: {r}")

            match r.restraint_type:
                case RestraintType.Covalent:
                    restraint_type = "covalent"

                case RestraintType.Pocket:
                    restraint_type = "pocket"
                case RestraintType.Contact:
                    restraint_type = "contact"
                case _:
                    raise ValueError(f"Unsupported restraint type: {r.restraint_type}")
            entries.append(
                (
                    chain1_id,
                    ...,
                    chain2_id,
                    ...,
                    restraint_type,
                    1.0,
                    0.0,
                    r.max_distance,
                    r.description,
                )
            )
        df = pd.DataFrame(entries).reset_index()
        df.columns = restraint_cols
        df.to_csv(self.restraints, index=False)

    def generate_chai_msa(self, msa_dir: str | Path):
        """Generate Chai-1 style MSA files."""
        ...

    def generate_chai_templates(self):
        """Generate Chai-1 style templates."""
        ...

    @staticmethod
    def _build_fasta_entry(
        seq_type: str, chain_id: str | list[str], sequence: str
    ) -> list[str]:
        """Build FASTA entry string."""
        if isinstance(chain_id, list):
            return [f">{seq_type}|{c}\n{sequence}" for c in chain_id]
        else:
            return [f">{seq_type}|{chain_id}\n{sequence}"]

    @staticmethod
    def _add_modifications_to_sequence(
        sequence: str, modifications: list[SequenceModification]
    ) -> str:
        """Add modifications to sequence string.

        When modified residues already have a CCD code, we should replace the residue
        with the CCD code in parentheses (no need to specify them as bonds).
        """
        if not modifications:
            return sequence
        seq_list = list(sequence)
        for mod in modifications:
            seq_list[mod.position - 1] = f"({mod.ccd})"
        return "".join(seq_list)

    @staticmethod
    def ccd_to_smiles(ccd: str, ccd_lib_dir: str | Path) -> str:
        """Convert CCD ligand to SMILES string."""
        import pickle

        from rdkit import Chem

        ccd_file = Path(ccd_lib_dir).expanduser().resolve() / f"{ccd}.pkl"
        with open(ccd_file, "rb") as f:
            mol = pickle.load(f)  # noqa: S301

        return Chem.MolToSmiles(mol)

    @staticmethod
    def _build_res_idx(restraint_type: str, chain_type: str, atom: Atom) -> str:
        """Build residue index string for Chai-1 restraints.

        Covalent bond: N436@N (residue), @C1 (ligand)
        Contact: R84, G7
        """
        if chain_type in {"ligand", "glycan"}:
            return f"@{atom.atom_name}"
        else:
            return f"{atom.residue_name}{atom.residue_idx}@{atom.atom_name}"


def run_chai(abcfold_conf_file: str | Path, output_dir: str | Path) -> bool:
    """Run Chai-1 using the ABCFold config file.

    Returns:
        Bool: True if the Chai-1 run was successful, False otherwise.

    """
    input_conf_file = Path(abcfold_conf_file).expanduser().resolve()
    conf = load_abcfold_config(input_conf_file)
    workdir = Path(output_dir).expanduser().resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    ...


# TODO: modify this to directly call chai_lab.chai1.run_inference
# def generate_chai_command(
#     input_fasta: str | Path,
#     msa_dir: str | Path,
#     input_constraints: str | Path,
#     output_dir: str | Path,
#     number_of_models: int = 5,
#     num_recycles: int = 10,
#     seed: int = 42,
#     use_templates_server: bool = False,
#     template_hits_path: Path | None = None,
# ) -> list:
#     """Generate the Chai-1 command

#     Args:
#         input_fasta (Union[str, Path]): Path to the input fasta file
#         msa_dir (Union[str, Path]): Path to the MSA directory
#         input_constraints (Union[str, Path]): Path to the input constraints file
#         output_dir (Union[str, Path]): Path to the output directory
#         number_of_models (int): Number of models to generate
#         num_recycles (int): Number of trunk recycles
#         seed (int): Seed for the random number generator
#         use_templates_server (bool): If True, use templates from the server
#         template_hits_path (Path): Path to the template hits m8 file

#     Returns:
#         list: The Chai-1 command

#     """
#     chai_exe = Path(__file__).parent / "chai.py"
#     cmd = ["python", str(chai_exe), "fold", str(input_fasta)]

#     if Path(msa_dir).exists():
#         cmd += ["--msa-directory", str(msa_dir)]
#     if Path(input_constraints).exists():
#         cmd += ["--constraint-path", str(input_constraints)]

#     cmd += ["--num-diffn-samples", str(number_of_models)]
#     cmd += ["--num-trunk-recycles", str(num_recycles)]
#     cmd += ["--seed", str(seed)]

#     assert not (use_templates_server and template_hits_path), (
#         "Cannot specify both templates server and path"
#     )

#     if shutil.which("kalign") is None and (use_templates_server or template_hits_path):
#         logger.warning(
#             "kalign not found, skipping template search kalign is required. \
# Please install kalign to use templates with Chai-1."
#         )
#     else:
#         if use_templates_server:
#             cmd += ["--use-templates-server"]
#         if template_hits_path:
#             cmd += ["--template-hits-path", str(template_hits_path)]

#     cmd += [str(output_dir)]

#     return cmd
