"""Entrypoint for ABCFold CLI."""

import typer

import abcfold.cli.prepare as cli_prepare

app = typer.Typer()

app.add_typer(cli_prepare.app, name="prepare", help="Prepare input files for folding.")


if __name__ == "__main__":
    app()
