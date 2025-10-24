"""Run Boltz using ABCFold config file."""

import logging
import subprocess as sp
import time
from datetime import UTC, datetime
from pathlib import Path

from tqdm import tqdm

from abcfold.schema import ABCFoldConfig

logger = logging.getLogger(__file__)


def run_boltz(
    abcfold_conf: ABCFoldConfig,
    output_dir: str | Path,
    boltz_yaml_file: str | Path,
    run_id: str,
) -> bool:
    """Entrypoint for running Boltz for structure prediction.

    Returns:
        True if the Boltz run was successful, False otherwise.

    """
    workdir = Path(output_dir).expanduser().resolve()
    log_dir = workdir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    boltz_yaml_file = Path(boltz_yaml_file).expanduser().resolve()
    if not boltz_yaml_file.exists():
        raise FileNotFoundError(f"Boltz config file not found: {boltz_yaml_file}")

    for seed in tqdm(abcfold_conf.seeds, desc=f"Boltz run {run_id}"):
        logger.info(f"Boltz run {run_id} using seed {seed}")
        cmd = generate_boltz_command(
            boltz_yaml_file,
            workdir / f"seed_{seed}",
            num_trunk_recycles=abcfold_conf.num_trunk_recycles,
            num_diffn_timesteps=abcfold_conf.num_diffn_timesteps,
            num_diffn_samples=abcfold_conf.num_diffn_samples,
            seed=seed,
            additional_args=abcfold_conf.boltz_additional_cli_args,
        )

        log_path = log_dir / f"{run_id}_boltz_seed{seed}.log"
        with (
            sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.STDOUT, encoding="utf-8") as p,  #  noqa: S603
            open(log_path, "w") as log_file,
        ):
            now = time.time()
            log_file.write(f"Time: {str(datetime.now(UTC))}\n")
            log_file.write(f"Running command: {' '.join(cmd)}\n\n")
            logger.info(f"Saving logs to {log_path}")

            logger.debug(f"Running command: {' '.join(cmd)}")
            stdout = ""
            while (buffered_output := p.stdout.readline()) != "" or p.poll() is None:
                stdout += buffered_output
                log_file.write(buffered_output)
                log_file.flush()

            log_file.write(f"\nFinished at: {str(datetime.now(UTC))}\n")
            log_file.write(f"Elapsed time: {time.time() - now:.2f} seconds\n")

            if p.returncode != 0:
                logger.error(f"Boltz run failed. Error log is in {log_path}")

                return False

            elif "WARNING: ran out of memory" in stdout:
                logger.error("Boltz ran out of memory")
                return False

    logger.info(f"Boltz run complete: {run_id}")
    logger.info(f"Output files are in {workdir}")
    return True


def generate_boltz_command(
    input_yaml: str | Path,
    output_dir: str | Path,
    num_trunk_recycles: int,
    num_diffn_timesteps: int,
    num_diffn_samples: int,
    additional_args: list[str] | None = None,
    seed: int = 42,
) -> list:
    """Generate the Boltz command."""
    out_path = Path(output_dir).expanduser().resolve()
    out_path.mkdir(parents=True, exist_ok=True)
    cmd = [
        "boltz",
        "predict",
        str(input_yaml),
        "--out_dir",
        str(out_path),
        "--recycling_steps",
        str(num_trunk_recycles),
        "--sampling_steps",
        str(num_diffn_timesteps),
        "--diffusion_samples",
        str(num_diffn_samples),
        "--seed",
        str(seed),
    ]
    if additional_args is not None:
        cmd.extend(additional_args)
    return cmd
