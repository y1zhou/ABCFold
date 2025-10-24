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
    force: Annotated[
        bool,
        typer.Option(
            "--force",
            "-f",
            help="Whether to overwrite existing MSA and template files.",
            flag_value=True,
        ),
    ] = False,
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
    """Query MSA and template files from ColabFold for all protein sequences.

    The generated files can be used by Chai directly.
    """
    conf_path = conf_file.expanduser().resolve()
    conf = load_abcfold_config(conf_path)
    out_path = out_dir.expanduser().resolve()
    if force:
        import shutil

        shutil.rmtree(out_path, ignore_errors=True)
    out_path.mkdir(parents=True, exist_ok=True)

    search_chains = set(chains.split(",")) if chains is not None else None
    conf = add_msa_to_config(
        conf,
        out_path / "msa",
        search_chains,
        search_templates,
        fetch_templates,
        template_cache_dir,
    )
    print(f"MSA files generated in: {out_path / 'msa'}")
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
    ccd_lib_dir: Annotated[
        Path | None,
        typer.Option(
            help="Path to the Boltz CCD library directory. Can omit if no CCD ligands or modifications were used.",
        ),
    ] = None,
):
    """Prepare FASTA and restraints files for Chai."""
    from abcfold.chai1.run_chai1_abcfold import ChaiConfig

    conf_path = conf_file.expanduser().resolve()
    conf = load_abcfold_config(conf_path)
    out_path = out_dir.expanduser().resolve()
    out_path.mkdir(parents=True, exist_ok=True)

    run_id = conf_path.stem
    chai_conf = ChaiConfig(conf, out_path, run_id, ccd_lib_dir)
    chai_conf.generate_chai_inputs()

    config_yaml_file = out_path / f"{run_id}.yaml"
    chai_conf.dump_chai_config(config_yaml_file)

    print(f"Chai FASTA file written to: {chai_conf.fasta}")
    if chai_conf.restraints is not None:
        print(f"Chai restraints file written to: {chai_conf.restraints}")
    print(f"Chai config YAML (with chain mappings) written to: {config_yaml_file}")


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
    max_num_templates_per_chain: Annotated[
        int,
        typer.Option(
            help="Maximum number of templates to use per chain. The default is 4, which matches Chai's limit. A higher number may lead to excessive GPU memory usage.",
        ),
    ] = 4,
):
    """Prepare YAML and MSA files for Boltz."""
    from abcfold.boltz.schema import abcfold_to_boltz

    conf_path = conf_file.expanduser().resolve()
    conf = load_abcfold_config(conf_path)
    out_path = out_dir.expanduser().resolve()
    out_path.mkdir(parents=True, exist_ok=True)

    boltz_conf = abcfold_to_boltz(
        conf, out_path / "boltz_msa", max_num_templates_per_chain
    )
    run_id = conf_path.stem
    boltz_yaml_file = out_path / f"{run_id}.yaml"
    write_config(boltz_conf, boltz_yaml_file)
    print(f"Boltz config written to: {boltz_yaml_file}")


if __name__ == "__main__":
    app()
