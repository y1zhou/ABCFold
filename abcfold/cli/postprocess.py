"""Run the post-processing steps for ABCFold."""

import logging
import shutil
from pathlib import Path
from typing import Annotated, Any

import typer
from tqdm import tqdm

from abcfold.output.boltz import BoltzOutput
from abcfold.output.chai import ChaiOutput

logger = logging.getLogger(__name__)

app = typer.Typer()


def collect_models(
    af3_results_dir: Path | None,
    boltz_results_dir: Path | None,
    chai_results_dir: Path | None,
) -> tuple[None, BoltzOutput | None, ChaiOutput | None]:
    """Collect generated structure models and corresponding score files."""
    ao, bo, co = None, None, None
    if af3_results_dir is not None:
        raise NotImplementedError(
            "AF3 is very low priority. Please use `abcfold` instead."
        )

    if boltz_results_dir is not None:
        logger.info("Collecting Boltz results...")

        boltz_out_dirs = list(boltz_results_dir.glob("boltz_results_seed-*"))
        bo = PatchedBoltzOutput(boltz_out_dirs)

    if chai_results_dir is not None:
        logger.info("Collecting Chai results...")

        chai_out_dirs = list(chai_results_dir.glob("chai_seed-*"))
        co = PatchedChaiOutput(chai_out_dirs)

    return ao, bo, co


def cif_to_pdb_models(
    models: list[BoltzOutput | ChaiOutput],
    pdb_out_dir: str | Path,
    superpose_chains: str | None,
) -> list[BoltzOutput | ChaiOutput]:
    """Convert structure models to separate folder in PDB format.

    This is for maximum compatibility with pae-viewer, as it only has limited support
    for PDBx/mmCIF files.
    """
    from abcfold.output.file_handlers import CifFile, superpose_models

    pdb_out_path = Path(pdb_out_dir)

    pdb_model_paths: list[Path] = []
    logger.info("Converting output model files to PDB...")
    for m in models:
        for cif_models in m.cif_files.values():
            for cif_model in cif_models:
                cif_model: CifFile
                output_pdb_name = f"{cif_model.name}.pdb"
                output_model_path = pdb_out_path / output_pdb_name

                if not output_model_path.exists():
                    cif_model.model.write_pdb(str(output_model_path))

                pdb_model_paths.append(output_model_path)
                cif_model.pathway = output_model_path

    if len(pdb_model_paths) > 1:
        logger.info("Superpositioning output models...")
        superpose_models(pdb_model_paths, superpose_chains)

    return models


def build_output_pages(
    models: list[BoltzOutput | ChaiOutput], plot_paths: dict[str, str], out_path: Path
) -> tuple[list[str], dict[str, Any]]:
    """Combine models from different methods for post-processing."""
    from abcfold.html.html_utils import (
        get_all_cif_files,
        get_model_data,
        get_model_sequence_data,
    )
    from abcfold.output.file_handlers import superpose_models
    from abcfold.output.utils import get_gap_indicies, insert_none_by_minus_one

    cif_models = [
        cif_file
        for cif_list in get_all_cif_files(models).values()
        for cif_file in cif_list
    ]
    indicies = get_gap_indicies(*cif_models)
    index_counter = 0

    programs_run = []
    combined_models = []
    for m_output in models:
        if isinstance(m_output, BoltzOutput):
            program = "Boltz"
        elif isinstance(m_output, ChaiOutput):
            program = "Chai-1"
        else:
            # TODO: Add AlphaFold 3 post-processing when needed
            raise TypeError("Unknown model output type")

        programs_run.append(program)
        logger.info(f"Post-processing {program} models...")
        for seed in m_output.output.keys():
            for idx in m_output.output[seed].keys():
                model = m_output.output[seed][idx]["cif"]
                model.check_clashes()
                score_file = m_output.output[seed][idx]["scores"]
                plddt = model.residue_plddts
                if len(indicies) > 0:
                    plddt = insert_none_by_minus_one(indicies[index_counter], plddt)
                index_counter += 1
                model_data = get_model_data(
                    model, plot_paths, program, plddt, score_file, out_path
                )
                combined_models.append(model_data)

    logger.info("Preparing output score files...")
    sequence_data = get_model_sequence_data(cif_models)
    sequence = ""
    for key in sequence_data.keys():
        sequence += sequence_data[key]
    chain_data = {}
    ref = 0
    for key in sequence_data.keys():
        chain_data["Chain " + key] = (ref, len(sequence_data[key]) + ref - 1)
        ref += len(sequence_data[key])
    results_dict = {
        "sequence": sequence,
        "models": combined_models,
        "plotly_path": Path(plot_paths["plddt"]).relative_to(out_path).as_posix(),
        "chain_data": chain_data,
    }
    return programs_run, results_dict


@app.command(name="collect")
def postprocess(
    out_dir: Annotated[
        Path,
        typer.Option(
            ..., "--out-dir", "-o", help="Output directory for post-processed results."
        ),
    ],
    af3_results_dir: Annotated[
        Path | None,
        typer.Option(
            "--af3-results-dir",
            "-a",
            help="Directory containing AlphaFold 3 results to be post-processed.",
        ),
    ] = None,
    boltz_results_dir: Annotated[
        Path | None,
        typer.Option(
            "--boltz-results-dir",
            "-b",
            help="Directory containing Boltz results to be post-processed.",
        ),
    ] = None,
    chai_results_dir: Annotated[
        Path | None,
        typer.Option(
            "--chai-results-dir",
            "-c",
            help="Directory containing Chai results to be post-processed.",
        ),
    ] = None,
    superpose_chains: Annotated[
        str | None,
        typer.Option(
            help="Comma-separated chain IDs to use for superposition (e.g., 'A,B')."
            " Defaults to all chains.",
        ),
    ] = None,
):
    """Run post-processing steps for ABCFold.

    Adapted from abcfold.abcfold:run.
    """
    # Collect output models from different methods
    ao, bo, co = collect_models(af3_results_dir, boltz_results_dir, chai_results_dir)
    found_models = [x for x in (ao, bo, co) if x is not None]
    if not found_models:
        logger.warning("No output models found for further processing.")
        return

    out_path = Path(out_dir).expanduser().resolve()
    (out_path / "output_models").mkdir(parents=True, exist_ok=True)
    found_models = cif_to_pdb_models(
        found_models, out_path / "output_models", superpose_chains
    )

    # Compile data to make output pages
    from abcfold.abcfold import PLOTS_DIR
    from abcfold.html.html_utils import plots

    logger.info("Generating plots...")
    plot_dict = plots(found_models, out_path / PLOTS_DIR)
    programs_run, results_dict = build_output_pages(found_models, plot_dict, out_path)

    # Create the index page
    import orjson

    from abcfold.abcfold import HTML_DIR, HTML_TEMPLATE
    from abcfold.html.html_utils import render_template

    logger.info("Generating output HTML page...")
    results_json = orjson.dumps(results_dict).decode()
    if not (out_path / "scores.json").exists():
        with open(out_path / "scores.json", "w") as f:
            f.write(results_json)
    if not out_path.joinpath(".feature_viewer").exists():
        shutil.copytree(HTML_DIR, out_path / ".feature_viewer")

    programs = (
        f"Structure predictions for: {', '.join(programs_run[:-1])}, and {programs_run[-1]}"
        if len(programs_run) > 1
        else f"Structure predictions for: {programs_run[0]}"
    )

    html_out_path = out_path.joinpath("index.html").resolve()
    render_template(
        HTML_TEMPLATE,
        html_out_path,
        # kwargs appear as variables in the template
        abcfold_html_dir=".feature_viewer",
        programs=programs,
        results_json=results_json,
        version=0.1,
    )
    logger.info(f"Output page written to {html_out_path}")


@app.command(name="serve")
def serve(
    result_dir: Annotated[
        Path,
        typer.Argument(help="Directory containing post-processed results."),
    ],
    host: Annotated[
        str,
        typer.Option("--host", "-h", help="Host address to bind the server to."),
    ] = "localhost",
    port: Annotated[
        int,
        typer.Option("--port", "-p", help="Port number to bind the server to."),
    ] = 8000,
):
    """Serve the output HTML page with a local web server."""
    import os
    import socketserver
    import sys
    import webbrowser

    from abcfold.html.html_utils import (
        NoCacheHTTPRequestHandler,
        output_open_html_script,
    )

    result_path = Path(result_dir).expanduser().resolve()
    if not result_path.exists():
        logger.error(f"Result directory {result_path} does not exist.")
        return
    # Change to the output directory to run the server
    os.chdir(result_path)

    # Make a script to open the output HTML file in the default web browser
    output_open_html_script("open_output.py", port=port)

    try:
        # Start the server
        with socketserver.TCPServer(("", port), NoCacheHTTPRequestHandler) as httpd:
            logger.info(f"Serving at http://{host}:{port}/index.html")
            logger.info("Press Ctrl+C to stop the server")
            # Open the main HTML page in the default web browser
            webbrowser.open(f"http://{host}:{port}/index.html")
            # Keep the server running
            httpd.serve_forever()
    except KeyboardInterrupt:
        logger.info("Server stopped")
        httpd.server_close()
        sys.exit(0)


class PatchedBoltzOutput(BoltzOutput):
    """Patched BoltzOutput to bypass moving files around.

    We already have the files in place from ABCFold.
    """

    def __init__(self, boltz_output_dirs: list[str | Path]):
        """Process the output of a Boltz run.

        Ref: `abcfold.abcfold.output.boltz:BoltzOutput`
        """
        self.output_dirs = [Path(x) for x in boltz_output_dirs]

        # Placeholder attributes to satisfy the parent class
        self.input_params = None
        self.yaml_input_obj = None

        self.output = self.process_boltz_output()
        self.seeds = list(self.output.keys())
        self.pae_files = {
            seed: [value["pae"] for value in self.output[seed].values()]
            for seed in self.seeds
        }
        self.cif_files = {
            seed: [value["cif"] for value in self.output[seed].values()]
            for seed in self.seeds
        }
        self.plddt_files = {
            seed: [value["plddt"] for value in self.output[seed].values()]
            for seed in self.seeds
        }
        self.pde_files = {
            seed: [value["pde"] for value in self.output[seed].values()]
            for seed in self.seeds
        }
        self.scores_files = {
            seed: [value["json"] for value in self.output[seed].values()]
            for seed in self.seeds
        }
        self.pae_to_af3()
        self.af3_pae_files = {
            seed: [value["af3_pae"] for value in self.output[seed].values()]
            for seed in self.seeds
        }

    def update_chain_labels(self, cif_file):
        """Boltz generated files have matching chain labels to the input.

        Called when parsing files in `self.process_boltz_output()`.
        """
        return cif_file

    def pae_to_af3(self):
        """Convert the PAE data from Boltz to the format used by Alphafold3.

        Avoid modifying files in-place; instead, write files with a '-af3.json' suffix.
        """
        from abcfold.output.file_handlers import ConfidenceJsonFile
        from abcfold.output.utils import Af3Pae

        new_pae_files = {}
        for seed in tqdm(self.seeds, desc="Converting Boltz PAE to AF3 format"):
            for pae_file, cif_file in zip(
                self.pae_files[seed], self.cif_files[seed], strict=True
            ):
                pae_path = pae_file.pathway
                out_name = pae_path.with_name(f"{pae_path.stem}-af3.json")
                if not out_name.exists():
                    pae = Af3Pae.from_boltz(pae_file.data, cif_file)
                    pae.to_file(out_name)

                if seed not in new_pae_files:
                    new_pae_files[seed] = []
                new_pae_files[seed].append(ConfidenceJsonFile(out_name))

        self.output = {
            seed: {
                i: {
                    "cif": cif_file,
                    "af3_pae": new_pae_files[seed][i],
                    "scores": self.output[seed][i]["json"],
                }
                for i, cif_file in enumerate(self.cif_files[seed])
            }
            for seed in self.seeds
        }


class PatchedChaiOutput(ChaiOutput):
    """Patched ChaiOutput to bypass get_input_fasta()."""

    def __init__(self, chai_output_dirs: list[str | Path]):
        """Avoid moving files around for Chai output processing."""
        self.output_dirs = [Path(x) for x in chai_output_dirs]

        # Placeholder attributes to satisfy the parent class
        self.input_params = None
        self.input_fasta = None

        self.output = self.process_chai_output()
        self.seeds = list(self.output.keys())

        self.pae_files = {
            seed: [
                value["pae"] for value in self.output[seed].values() if "pae" in value
            ]
            for seed in self.seeds
        }
        self.cif_files = {
            seed: [value["cif"] for value in self.output[seed].values()]
            for seed in self.seeds
        }
        self.scores_files = {
            seed: [value["scores"] for value in self.output[seed].values()]
            for seed in self.seeds
        }
        self.pae_to_af3()
        self.af3_pae_files = {
            seed: [value["af3_pae"] for value in self.output[seed].values()]
            for seed in self.seeds
        }

    def process_chai_output(self):
        """Process the output of a Chai run.

        Main difference from parent class:

        - Bypass calling `self.update_chain_labels()`.
        - Only glob for '*.model_idx*' files instead of all possible output files.
        """
        import numpy as np

        from abcfold.output.file_handlers import CifFile, FileTypes, NpyFile, NpzFile

        file_groups = {}

        for pathway in tqdm(self.output_dirs, desc="Processing Chai output"):
            seed = pathway.name.split("_")[-1]
            if seed not in file_groups:
                file_groups[seed] = {}

            # Handle PAE scores
            if (pae_file := pathway / "pae_scores.npy").exists():
                # Compatibility with previous implementation which dumped all PAE scores into one file
                # When such cases are detected, split them into per-model files and remove the old file
                pae_scores = np.load(pae_file)
                for i in range(pae_scores.shape[0]):
                    np.save(pathway / f"pae.model_idx_{i}.npy", pae_scores[i])
                pae_file.unlink()

            for output in pathway.rglob("*.model_idx_*"):
                number = output.stem.split(".model_idx_")[-1]
                if number.isdigit():
                    number = int(number)
                else:
                    number = -1

                file_type = output.suffix[1:]
                match file_type:
                    case FileTypes.NPZ.value:
                        file_ = NpzFile(str(output))

                    case FileTypes.CIF.value:
                        file_ = CifFile(str(output), self.input_params)
                        # file_ = self.update_chain_labels(file_)

                    case FileTypes.NPY.value:
                        file_ = NpyFile(str(output))
                    case _:
                        continue

                if number not in file_groups[seed]:
                    file_groups[seed][number] = [file_]
                else:
                    file_groups[seed][number].append(file_)

        seed_dict = {}
        for seed, models in tqdm(file_groups.items(), desc="Collecting Chai scores"):
            model_number_file_type_file = {}
            pae_file_data = None
            if -1 in models:
                for file_ in models[-1]:
                    if file_.pathway.name == "pae_scores.npy":
                        pae_file_data = file_
                        break

            for model_number, files in models.items():
                if model_number == -1:
                    continue
                intermediate_dict = {}
                for file_ in sorted(files, key=lambda x: x.suffix):
                    if file_.pathway.stem.startswith("scores.model"):
                        intermediate_dict["scores"] = file_
                    elif file_.pathway.stem.startswith("pred.model"):
                        file_.name = f"Chai-1_{seed}_{model_number}"
                        # Chai cif not recognised by pae-viewer, so we load and save
                        # TODO: pae-viewer requires the original .to_file?
                        # file_.pathway = (
                        #     file_.pathway.parent / f"stripped-{file_.pathway.name}"
                        # )
                        # file_.to_file(file_.pathway)
                        intermediate_dict["cif"] = file_
                    elif file_.pathway.stem.startswith("pae.model"):
                        intermediate_dict["pae"] = file_

                model_number_file_type_file[model_number] = intermediate_dict

            model_number_file_type_file = {
                model_number: model_number_file_type_file[model_number]
                for model_number in sorted(model_number_file_type_file)
            }
            seed_dict[seed] = model_number_file_type_file

        return seed_dict

    def pae_to_af3(self):
        """Convert the Chai-1 PAE data to the format expected by AlphaFold3.

        Avoid modifying files in-place; instead, write files with a '-af3.json' suffix.
        """
        from abcfold.output.file_handlers import ConfidenceJsonFile
        from abcfold.output.utils import Af3Pae

        new_pae_files = {}
        for seed in tqdm(self.seeds, desc="Converting Chai PAE to AF3 format"):
            for i, (pae_file, cif_file) in enumerate(
                zip(self.pae_files[seed], self.cif_files[seed], strict=True)
            ):
                pae_path = pae_file.pathway
                out_name = pae_path.with_name(f"{pae_path.stem}-af3.json")
                if not out_name.exists():
                    pae = Af3Pae.from_chai1(pae_file.data, cif_file)
                    pae.to_file(out_name)

                if seed not in new_pae_files:
                    new_pae_files[seed] = []
                new_pae_files[seed].append(ConfidenceJsonFile(out_name))

        self.output = {
            seed: {
                i: {
                    "cif": cif_file,
                    "af3_pae": new_pae_files[seed][i],
                    "scores": self.output[seed][i]["scores"],
                }
                for i, cif_file in enumerate(self.cif_files[seed])
            }
            for seed in self.seeds
        }
