import logging
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Union

from abcfold.chai1.af3_to_chai import ChaiFasta

logger = logging.getLogger("logger")


def run_chai(
    input_json: Union[str, Path],
    output_dir: Union[str, Path],
    save_input: bool = False,
    test: bool = False,
    number_of_models: int = 5,
    num_recycles: int = 10,
    use_templates_server: bool = False,
    template_hits_path: Path | None = None,
) -> bool:
    """
    Run Chai-1 using the input JSON file

    Args:
        input_json (Union[str, Path]): Path to the input JSON file
        output_dir (Union[str, Path]): Path to the output directory
        save_input (bool): If True, save the input fasta file and MSA to the output
        directory
        test (bool): If True, run the test command
        number_of_models (int): Number of models to generate
        num_recycles (int): Number of trunk recycles
        use_templates_server (bool): If True, use templates from the server
        template_hits_path (Path): Path to the template hits m8 file

    Returns:
        Bool: True if the Chai-1 run was successful, False otherwise

    """
    input_json = Path(input_json)
    output_dir = Path(output_dir)

    # logger.debug("Checking if Chai-1 is installed")
    # check_chai1()

    with tempfile.TemporaryDirectory() as temp_dir:
        working_dir = Path(temp_dir)
        if save_input:
            logger.info("Saving input fasta file and msa to the output directory")
            working_dir = output_dir
            working_dir.mkdir(parents=True, exist_ok=True)

        chai_fasta = ChaiFasta(working_dir)
        chai_fasta.json_to_fasta(input_json)

        out_fasta = chai_fasta.fasta
        msa_dir = chai_fasta.working_dir
        out_constraints = chai_fasta.constraints

        for seed in chai_fasta.seeds:
            chai_output_dir = output_dir / f"chai_output_seed-{seed}"

            logger.info(f"Running Chai-1 using seed {seed}")
            cmd = (
                generate_chai_command(
                    out_fasta,
                    msa_dir,
                    out_constraints,
                    chai_output_dir,
                    number_of_models,
                    num_recycles=num_recycles,
                    seed=seed,
                    use_templates_server=use_templates_server,
                    template_hits_path=template_hits_path,
                )
                if not test
                else generate_chai_test_command()
            )
            logger.debug("Running command: %s", " ".join(cmd))
            with subprocess.Popen(
                cmd,
                stdout=sys.stdout,
                stderr=subprocess.PIPE,
            ) as proc:
                _, stderr = proc.communicate()
                if proc.returncode != 0:
                    if proc.stderr:
                        if chai_output_dir.exists():
                            output_err_file = chai_output_dir / "chai_error.log"
                        else:
                            output_err_file = chai_output_dir.parent / "chai_error.log"
                        with open(output_err_file, "w") as f:
                            f.write(stderr.decode())
                        logger.error(
                            "Chai-1 run failed. Error log is in %s", output_err_file
                        )
                    else:
                        logger.error("Chai-1 run failed")
                    return False

        logger.info("Chai-1 run complete")
        return True


def generate_chai_command(
    input_fasta: Union[str, Path],
    msa_dir: Union[str, Path],
    input_constraints: Union[str, Path],
    output_dir: Union[str, Path],
    number_of_models: int = 5,
    num_recycles: int = 10,
    seed: int = 42,
    use_templates_server: bool = False,
    template_hits_path: Path | None = None,
) -> list:
    """
    Generate the Chai-1 command

    Args:
        input_fasta (Union[str, Path]): Path to the input fasta file
        msa_dir (Union[str, Path]): Path to the MSA directory
        input_constraints (Union[str, Path]): Path to the input constraints file
        output_dir (Union[str, Path]): Path to the output directory
        number_of_models (int): Number of models to generate
        num_recycles (int): Number of trunk recycles
        seed (int): Seed for the random number generator
        use_templates_server (bool): If True, use templates from the server
        template_hits_path (Path): Path to the template hits m8 file

    Returns:
        list: The Chai-1 command

    """

    chai_exe = Path(__file__).parent / "chai.py"
    cmd = ["python", str(chai_exe), "fold", str(input_fasta)]

    if Path(msa_dir).exists():
        cmd += ["--msa-directory", str(msa_dir)]
    if Path(input_constraints).exists():
        cmd += ["--constraint-path", str(input_constraints)]

    cmd += ["--num-diffn-samples", str(number_of_models)]
    cmd += ["--num-trunk-recycles", str(num_recycles)]
    cmd += ["--seed", str(seed)]

    assert not (use_templates_server and template_hits_path), (
        "Cannot specify both templates server and path"
    )

    if shutil.which("kalign") is None and (use_templates_server or template_hits_path):
        logger.warning(
            "kalign not found, skipping template search kalign is required. \
Please install kalign to use templates with Chai-1."
        )
    else:
        if use_templates_server:
            cmd += ["--use-templates-server"]
        if template_hits_path:
            cmd += ["--template-hits-path", str(template_hits_path)]

    cmd += [str(output_dir)]

    return cmd


def generate_chai_test_command() -> list:
    """
    Generate the Chai-1 test command

    Args:
        None

    Returns:
        list: The Chai-1 test command
    """
    chai_exe = Path(__file__).parent / "chai.py"
    return [
        "python",
        chai_exe,
        "fold",
        "--help",
    ]
