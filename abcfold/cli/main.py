"""Entrypoint for ABCFold CLI."""

from pathlib import Path
from typing import Annotated

import typer
from rich import print_json

import abcfold.cli.prepare as cli_prepare
from abcfold.schema import load_abcfold_config

app = typer.Typer()

app.add_typer(cli_prepare.app, name="prepare", help="Prepare input files for folding.")


@app.command(name="validate")
def validate_config(
    config_file: Annotated[
        Path, typer.Argument(..., help="Path to the ABCFold configuration file.")
    ],
):
    """Validate an ABCFold configuration file."""
    conf = load_abcfold_config(config_file)
    print_json(
        data=conf.model_dump(
            exclude_unset=True, exclude_none=True, exclude_computed_fields=True
        )
    )


if __name__ == "__main__":
    app()
