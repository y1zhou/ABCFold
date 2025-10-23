"""Generate input files for AF3, Chai, and Boltz."""

from pathlib import Path
from typing import Annotated

import typer

from abcfold.boltz.schema import abcfold_to_boltz
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
    chains: Annotated[
        str | None,
        typer.Option(
            help="Chains to process, e.g., A,B,C. Defaults to all protein chains."
        ),
    ] = None,
    search_templates: Annotated[
        bool, typer.Option(help="Whether to search for templates.", flag_value=True)
    ] = True,
    fetch_templates: Annotated[
        bool,
        typer.Option(help="Whether to fetch templates from RCSB.", flag_value=True),
    ] = True,
    template_cache_dir: Annotated[
        Path | None,
        typer.Option(
            help="Directory to cache fetched templates. Defaults to ~/.cache/rcsb/.",
        ),
    ] = None,
):
    """Query MSA and template files from ColabFold for all protein sequences."""
    conf_path = conf_file.expanduser().resolve()
    conf = load_abcfold_config(conf_path)
    out_path = out_dir.expanduser().resolve()
    out_path.mkdir(parents=True, exist_ok=True)

    search_chains = set(chains.split(",")) if chains is not None else None
    conf = add_msa_to_config(
        conf,
        out_path,
        search_chains,
        search_templates,
        fetch_templates,
        template_cache_dir,
    )
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
    raise NotImplementedError("AF3 is very low priority. Please use `abcfold` instead.")


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
    conf_path = conf_file.expanduser().resolve()
    conf = load_abcfold_config(conf_path)
    out_path = out_dir.expanduser().resolve()
    out_path.mkdir(parents=True, exist_ok=True)

    boltz_conf = abcfold_to_boltz(conf, out_path / "boltz_msa")
    run_id = conf_path.stem
    boltz_yaml_file = out_path / f"{run_id}.yaml"
    write_config(boltz_conf, boltz_yaml_file)
    print(f"Boltz config written to: {boltz_yaml_file}")


if __name__ == "__main__":
    app()
