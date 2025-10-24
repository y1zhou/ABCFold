"""Run Chai-1 using ABCFold configuration file."""

import logging
import time
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd
import yaml
from tqdm import tqdm

from abcfold.schema import (
    ABCFoldConfig,
    Atom,
    Glycan,
    Ligand,
    Polymer,
    ProteinSeq,
    RestraintType,
    SequenceModification,
)

logger = logging.getLogger("logger")


class ChaiConfig:
    """Builder for Chai-1 config files from ABCFold config."""

    def __init__(
        self,
        conf: ABCFoldConfig,
        work_dir: str | Path,
        run_id: str,
        ccd_lib_dir: str | Path | None = None,
    ):
        """Initialize ChaiConfig builder.

        Args:
            conf (ABCFoldConfig): ABCFold configuration object.
            work_dir (str | Path): Output directory for Chai-1 inputs.
            run_id (str): Unique identifier for the run.
            ccd_lib_dir: Boltz cache directory for CCD ligands.

        """
        self.conf = conf
        self.work_dir = Path(work_dir).expanduser().resolve()
        self.run_id = run_id
        self.ccd_lib_dir = (
            Path(ccd_lib_dir).expanduser().resolve()
            if ccd_lib_dir is not None
            else None
        )

        self.work_dir.mkdir(parents=True, exist_ok=True)

    def generate_chai_inputs(self):
        """Generate all Chai-1 input files."""
        # File paths
        self.fasta = self.work_dir / f"{self.run_id}.fasta"
        self.restraints: Path | None = (
            self.work_dir / f"{self.run_id}.restraints"
            if self.conf.restraints
            else None
        )
        self.msa: Path | None = self.work_dir / "chai_msa"
        # TODO: handle templates

        # Other metadata
        self.chain_type: dict[str, str] = {}
        self.fasta_entries = self.generate_chai_fasta(self.ccd_lib_dir)
        self.chain_id = self._build_chain_id_mapping()
        self.restraints_df = self.generate_chai_restraints()
        self.generate_chai_msa()

    def generate_chai_fasta(self, ccd_lib_dir: str | Path | None = None) -> list[str]:
        """Build Chai-1 style FASTA from ABCFold config.

        `ccd_lib_dir` is needed if their are CCD ligands in the input.
        Chai only supports SMILES ligands.

        Note that Chai-1 outputs re-assigns chain IDs based on the order of sequences
        in the FASTA file. So the chain IDs in the ABCFold config may not be preserved.

        Returns:
            List of FASTA entries. This method also modifies self.chain_type_map.

        """
        entries: list[str] = []
        for seq in self.conf.sequences:
            if isinstance(seq, Glycan):
                for chain_id in seq.id if isinstance(seq.id, list) else [seq.id]:
                    self.chain_type[chain_id] = "glycan"
                entries.extend(self._build_fasta_entry("glycan", seq.id, seq.chai_str))

            elif isinstance(seq, Ligand):
                for chain_id in seq.id if isinstance(seq.id, list) else [seq.id]:
                    self.chain_type[chain_id] = "ligand"
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
                for chain_id in seq.id if isinstance(seq.id, list) else [seq.id]:
                    self.chain_type[chain_id] = seq.seq_type.value
                entries.extend(
                    self._build_fasta_entry(
                        seq.seq_type.value,
                        seq.seq_hash,  # need to match query_id in m8 templates file
                        self._add_modifications_to_sequence(
                            seq.sequence, seq.modifications
                        ),
                    )
                )
            else:
                raise ValueError(f"Unsupported sequence type: {type(seq)}")

        with open(self.fasta, "w") as f:
            f.write("\n".join(entries) + "\n")

        return entries

    def generate_chai_restraints(self) -> pd.DataFrame | None:
        """Generate Chai-1 style restraints CSV from ABCFold config.

        Ref: https://github.com/chaidiscovery/chai-lab/issues/326#issuecomment-2862206725
        Ref: https://github.com/chaidiscovery/chai-lab/tree/main/examples/covalent_bonds
        Ref: https://github.com/chaidiscovery/chai-lab/tree/main/examples/restraints
        """
        if not self.conf.restraints:
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
            # Get Chai-assigned chain IDs
            chain1_id = self.chain_id[r.atom1.chain_id]
            chain2_id = self.chain_id[r.atom2.chain_id]
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

            res1_idx = self._build_res_idx(
                restraint_type,
                self.chain_type[r.atom1.chain_id],
                r.atom1,
            )
            res2_idx = self._build_res_idx(
                restraint_type,
                self.chain_type[r.atom2.chain_id],
                r.atom2,
            )

            # Special treatment for pocket restraints: swap to put hotspot in chain B
            if restraint_type == "pocket":
                if chain2_id == r.boltz_binder_chain:
                    chain1_id, chain2_id, res1_idx, res2_idx = (
                        chain2_id,
                        chain1_id,
                        None,
                        res1_idx,
                    )
                elif chain1_id == r.boltz_binder_chain:
                    res1_idx = None
                else:
                    raise ValueError(
                        f"Pocket restraint binder chain ID {r.boltz_binder_chain} "
                        f"not found in restraint atoms: {r}"
                    )
            entries.append(
                (
                    chain1_id,
                    res1_idx,
                    chain2_id,
                    res2_idx,
                    restraint_type,
                    1.0,
                    0.0,
                    r.max_distance,
                    r.description,
                )
            )
            # TODO: validate residue names == Chai-1 expectations
        df = pd.DataFrame(entries).reset_index()
        df.columns = restraint_cols
        df.to_csv(self.restraints, index=False)
        return df

    def generate_chai_msa(self):
        """Generate Chai-1 style MSA directory from ABCFold config."""
        protein_chains = [
            seq for seq in self.conf.sequences if isinstance(seq, ProteinSeq)
        ]

        # No need to move stuff around if there's only one shared MSA directory
        # This should be the case when the MSAs were prepared using ABCFold CLI
        msa_dirs = {c.msa_dir for c in protein_chains if c.msa_dir is not None}
        if not msa_dirs:
            self.msa = None
            return
        if len(msa_dirs) == 1:
            self.msa = Path(msa_dirs.pop()).expanduser().resolve()
            return

        # Link `.aligned.pqt` and `a3ms/*.a3m` files
        (self.msa / "a3ms").mkdir(parents=True, exist_ok=True)
        for msa_dir in msa_dirs:
            msa_dir_path = Path(msa_dir).expanduser().resolve()
            for msa_file in msa_dir_path.glob("*.aligned.pqt"):
                target_file = self.msa / msa_file.name
                if not target_file.exists():
                    target_file.symlink_to(msa_file)

            for a3m_file in (msa_dir_path / "a3ms").glob("*.a3m"):
                target_file = self.msa / "a3ms" / a3m_file.name
                if not target_file.exists():
                    target_file.symlink_to(a3m_file)

    def dump_chai_config(self, file_path: str | Path):
        """Dump Chai-1 config summary to a YAML file.

        The YAML file contains the following fields:

        - run_id: Unique identifier for the run.
        - fasta: Path to the Chai FASTA file.
        - restraints: Path to the restraints file (if any).
        - msa: Path to the MSA directory (if any).
        - chain_id_mapping: Mapping from original chain IDs to Chai-assigned chain IDs.
        """
        chai_conf = {
            "run_id": self.run_id,
            "fasta": str(self.fasta),
            "restraints": str(self.restraints) if self.restraints else None,
            "msa": str(self.msa) if self.msa else None,
            "chain_id_mapping": self.chain_id,
        }
        with open(file_path, "w") as f:
            yaml.safe_dump(chai_conf, f, default_flow_style=None)

    def _build_fasta_entry(
        self, seq_type: str, chain_id: str | list[str], sequence: str
    ) -> list[str]:
        """Build FASTA entry string."""
        if isinstance(chain_id, list):
            return [f">{seq_type}|{c}\n{sequence}" for c in chain_id]
        else:
            return [f">{seq_type}|{chain_id}\n{sequence}"]

    @staticmethod
    def _add_modifications_to_sequence(
        sequence: str, modifications: list[SequenceModification] | None
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

    def _build_chain_id_mapping(self) -> dict[str, str]:
        """Build mapping from original chain IDs to Chai-assigned chain IDs."""
        chai_chain_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        mapping = {}
        i = 0
        for seq in self.conf.sequences:
            for seq_id in seq.id if isinstance(seq.id, list) else [seq.id]:
                mapping[seq_id] = chai_chain_ids[i]
                i += 1
        return mapping

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
        # For ligands and glycans, only atom name is used
        if chain_type in {"ligand", "glycan"}:
            return f"@{atom.atom_name}"

        if restraint_type == "covalent":
            return f"{atom.residue_name}{atom.residue_idx}@{atom.atom_name}"
        else:  # contact or pocket
            return f"{atom.residue_name}{atom.residue_idx}"


def run_chai(
    # ABCFold configuration
    abcfold_conf: ABCFoldConfig,
    output_dir: str | Path,
    chai_yaml_file: str | Path,
    # Configuration for ESM, MSA, constraints, and templates
    template_hits_path: Path | None = None,
    template_cif_dir: Path | None = None,
    use_esm_embeddings: bool = True,
    # Parameters controlling how we do inference
    recycle_msa_subsample: int = 0,
    low_memory: bool = True,
) -> bool:
    """Entrypoint for running Chai-1 structure prediction.

    Args:
        abcfold_conf: ABCFoldConfig object.
        output_dir: Output directory for Chai-1 results.
        chai_yaml_file: Path to Chai-1 config YAML file generated by
            `abcfold.cli.prepare:prepare_chai`.
        template_hits_path: Path to template hits file (in m8 format).
        template_cif_dir: Directory containing template mmCIF files.
        use_esm_embeddings: Whether to use ESM embeddings.
        recycle_msa_subsample: If >0, subsample the input MSA.
        low_memory: Whether to run Chai-1 in low memory mode.

    Returns:
        True if the Chai-1 run was successful, False otherwise.

    """
    import os
    from contextlib import redirect_stderr, redirect_stdout

    workdir = Path(output_dir).expanduser().resolve()
    log_dir = workdir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    chai_yaml_file = Path(chai_yaml_file).expanduser().resolve()
    if not chai_yaml_file.exists():
        raise FileNotFoundError(f"Chai-1 config YAML file not found: {chai_yaml_file}")
    with open(chai_yaml_file) as f:
        chai_conf = yaml.safe_load(f)

    if template_cif_dir is not None:
        os.environ["CHAI_TEMPLATE_CIF_FOLDER"] = str(template_cif_dir)

    from chai_lab.chai1 import run_inference

    run_id = chai_conf["run_id"]
    for seed in tqdm(abcfold_conf.seeds, desc=f"Chai-1 run {run_id}"):
        logger.info(f"Chai-1 run {run_id} using seed {seed}")
        log_path = log_dir / f"{run_id}_boltz_seed{seed}.log"
        logger.info(f"Saving logs to {log_path}")
        with (
            open(log_path, "w") as log_file,
            redirect_stdout(log_file),
            redirect_stderr(log_file),
        ):
            now = time.time()
            log_file.write(f"Time: {str(datetime.now(UTC))}\n")

            run_dir = workdir / f"seed_{seed}"
            run_dir.mkdir(parents=True, exist_ok=True)
            success = run_inference(
                fasta_file=Path(chai_conf["fasta"]),
                output_dir=run_dir,
                msa_directory=Path(chai_conf["msa"]) if chai_conf["msa"] else None,
                constraint_path=Path(chai_conf["restraints"])
                if chai_conf["restraints"]
                else None,
                use_esm_embeddings=use_esm_embeddings,
                template_hits_path=template_hits_path,
                recycle_msa_subsample=recycle_msa_subsample,
                low_memory=low_memory,
                seed=seed,
            )
            if not success:
                logger.error(f"Chai-1 run failed using seed {seed}.")
                return False

            log_file.write(f"\nFinished at: {str(datetime.now(UTC))}\n")
            log_file.write(f"Elapsed time: {time.time() - now:.2f} seconds\n")

    logger.info(f"Chai run complete: {run_id}")
    logger.info(f"Output files are in {workdir}")
    return True
