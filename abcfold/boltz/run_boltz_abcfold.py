"""Run Boltz using ABCFold config file."""

import logging
import subprocess as sp
import time
from datetime import UTC, datetime
from pathlib import Path

from abcfold.boltz.schema import abcfold_to_boltz
from abcfold.schema import load_abcfold_config, write_config

logger = logging.getLogger("logger")


def run_boltz(abcfold_conf_file: str | Path, output_dir: str | Path) -> bool:
    """Run Boltz using the ABCFold config file."""
    input_conf_file = Path(abcfold_conf_file).expanduser().resolve()
    conf = load_abcfold_config(input_conf_file)
    workdir = Path(output_dir).expanduser().resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    boltz_conf = abcfold_to_boltz(conf, workdir / "boltz_msa")
    run_id = input_conf_file.stem
    boltz_yaml_file = workdir / f"{run_id}.yaml"
    write_config(boltz_conf, boltz_yaml_file)

    log_dir = workdir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    for seed in conf.seeds:
        logger.info(f"Running Boltz using seed: {seed}")
        cmd = generate_boltz_command(
            boltz_yaml_file,
            workdir / "output",
            num_trunk_recycles=conf.num_trunk_recycles,
            num_diffn_timesteps=conf.num_diffn_timesteps,
            num_diffn_samples=conf.num_diffn_samples,
            seed=seed,
            additional_args=conf.boltz_additional_cli_args,
        )

        log_path = log_dir / f"{run_id}_boltz_seed{seed}.log"
        with (
            sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.STDOUT, encoding="utf-8") as p,  #  noqa: S603
            open(log_path, "w") as log_file,
        ):
            now = time.time()
            log_file.write(f"Time: {str(datetime.now(UTC))}\n")
            log_file.write(f"Running command: {' '.join(cmd)}\n\n")

            logger.debug(f"Running command: {' '.join(cmd)}\nSaving log to {log_path}")
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

    logger.info("Boltz run complete")
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
    cmd = [
        "boltz",
        "predict",
        str(input_yaml),
        "--out_dir",
        str(output_dir),
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
