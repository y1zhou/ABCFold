"""Run Boltz using ABCFold config file."""

import logging
import subprocess as sp
import time
from datetime import UTC, datetime
from pathlib import Path

logger = logging.getLogger(__file__)


def run_boltz(
    output_dir: str | Path,
    boltz_yaml_file: str | Path,
    seed: int,
    num_trunk_recycles: int = 3,  # recycling_steps
    num_diffn_timesteps: int = 200,  # sampling_steps
    num_diffn_samples: int = 5,  # diffusion_samples
    boltz_additional_cli_args: list[str] | None = None,
) -> Path:
    """Entrypoint for running Boltz for structure prediction.

    Returns:
        Path to the Boltz output directory.

    """
    workdir = Path(output_dir).expanduser().resolve()
    log_dir = workdir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    boltz_yaml_file = Path(boltz_yaml_file).expanduser().resolve()
    if not boltz_yaml_file.exists():
        raise FileNotFoundError(f"Boltz config file not found: {boltz_yaml_file}")

    run_name = boltz_yaml_file.stem
    cmd = generate_boltz_command(
        boltz_yaml_file,
        workdir / f"boltz_seed_{seed}",
        num_trunk_recycles=num_trunk_recycles,
        num_diffn_timesteps=num_diffn_timesteps,
        num_diffn_samples=num_diffn_samples,
        seed=seed,
        additional_args=boltz_additional_cli_args,
    )

    # Skip if output already exists
    final_out_dir = workdir / f"boltz_results_seed-{seed}"
    if all(
        (
            final_out_dir / "predictions" / run_name / f"pae_{run_name}_model_{i}.npz"
        ).exists()
        for i in range(num_diffn_samples)
    ):
        return final_out_dir

    log_path = log_dir / f"boltz_seed-{seed}.log"
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
            raise sp.CalledProcessError(p.returncode, cmd)

        elif "WARNING: ran out of memory" in stdout:
            raise MemoryError("Boltz ran out of memory")

    # Move Boltz output files out of their subdirectories
    (workdir / f"boltz_seed_{seed}" / f"boltz_results_{run_name}").rename(final_out_dir)
    (workdir / f"boltz_seed_{seed}").rmdir()

    return final_out_dir


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
