"""Run ABCFold with the given configuration."""

from pathlib import Path
from typing import Annotated

import typer

from abcfold.boltz.run_boltz_abcfold import run_boltz
from abcfold.schema import load_abcfold_config

app = typer.Typer()


@app.command(name="af3")
def fold_af3(
    conf_file: Annotated[
        Path, typer.Argument(help="Path to the ABCFold configuration file.")
    ],
    out_dir: Annotated[
        Path,
        typer.Option(
            ...,
            "--out-dir",
            "-o",
            help="Output directory for prepared files.",
        ),
    ],
):
    """Build structure models with AF3."""
    raise NotImplementedError("AF3 is very low priority. Please use `abcfold` instead.")


@app.command(name="chai")
def fold_chai(
    conf_file: Annotated[
        Path, typer.Argument(help="Path to the ABCFold configuration file.")
    ],
    out_dir: Annotated[
        Path,
        typer.Option(
            ...,
            "--out-dir",
            "-o",
            help="Output directory for prepared files.",
        ),
    ],
):
    """Build structure models with Chai."""
    raise NotImplementedError


@app.command(name="boltz")
def fold_boltz(
    conf_file: Annotated[
        Path, typer.Argument(help="Path to the ABCFold configuration file.")
    ],
    boltz_yaml_file: Annotated[
        Path,
        typer.Option(
            ..., "--boltz-yaml-file", "-i", help="Path to the Boltz configuration file."
        ),
    ],
    out_dir: Annotated[
        Path,
        typer.Option(
            ...,
            "--out-dir",
            "-o",
            help="Output directory for Boltz.",
        ),
    ],
):
    """Build structure models with Boltz."""
    conf_path = conf_file.expanduser().resolve()
    conf = load_abcfold_config(conf_path)

    boltz_conf_path = boltz_yaml_file.expanduser().resolve()
    if not boltz_conf_path.exists():
        raise FileNotFoundError(f"Boltz config file not found: {boltz_conf_path}")

    out_path = out_dir.expanduser().resolve()
    out_path.mkdir(parents=True, exist_ok=True)
    run_boltz(conf, boltz_conf_path, out_path, boltz_conf_path.stem)


if __name__ == "__main__":
    app()
