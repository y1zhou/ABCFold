"""Generate input files for AF3, Chai, and Boltz."""

from pathlib import Path
from typing import Annotated

import typer

from abcfold.schema import add_msa_to_config, load_abcfold_config, write_config

app = typer.Typer()


@app.command(name="msa")
def search_msa(
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
    search_templates: Annotated[
        bool, typer.Option(help="Whether to search for templates.", flag_value=True)
    ] = True,
):
    """Query MSA and template files from ColabFold for all protein sequences."""
    conf_path = conf_file.expanduser().resolve()
    conf = load_abcfold_config(conf_path)
    out_path = out_dir.expanduser().resolve()
    out_path.mkdir(parents=True, exist_ok=True)

    conf = add_msa_to_config(conf, out_path / "msas", search_templates)
    new_conf_path = out_path / f"{conf_path.stem}.yaml"
    write_config(conf, new_conf_path)
    print(f"Updated config written to: {new_conf_path}")


@app.command(name="af3")
def prepare_af3(
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
    """Prepare JSON and MSA files for AF3."""
    raise NotImplementedError


@app.command(name="chai")
def prepare_chai(
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
    """Prepare FASTA, restraints, MSA, and templates for Chai."""
    raise NotImplementedError


@app.command(name="boltz")
def prepare_boltz(
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
    """Prepare YAML and MSA files for Boltz."""
    raise NotImplementedError


if __name__ == "__main__":
    app()
