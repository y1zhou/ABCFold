"""Run the post-processing steps for ABCFold."""

import logging
from pathlib import Path
from typing import Annotated

import typer
from tqdm import tqdm

from abcfold.output.boltz import BoltzOutput
from abcfold.output.chai import ChaiOutput

logger = logging.getLogger(__name__)

app = typer.Typer()


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
):
    """Run post-processing steps for ABCFold.

    Adapted from abcfold.abcfold:run.
    """
    # Collect output models from different methods
    found_models = []
    if af3_results_dir is not None:
        raise NotImplementedError(
            "AF3 is very low priority. Please use `abcfold` instead."
        )

    if boltz_results_dir is not None:
        logger.info("Post-processing Boltz results...")

        boltz_out_dirs = list(boltz_results_dir.glob("boltz_results_*"))
        bo = PatchedBoltzOutput(boltz_out_dirs)
        found_models.append(bo)

    if chai_results_dir is not None:
        logger.info("Post-processing Chai results...")

        chai_out_dirs = list(chai_results_dir.glob("chai_*_seed-*"))
        co = PatchedChaiOutput(chai_out_dirs)
        found_models.append(co)

    # Compile data to make output page
    if not found_models:
        logger.warning("No output models found for further processing.")
        return

    from abcfold.abcfold import PLOTS_DIR
    from abcfold.html.html_utils import get_all_cif_files, get_model_data, plots
    from abcfold.output.utils import get_gap_indicies, insert_none_by_minus_one

    out_path = Path(out_dir).expanduser().resolve()
    out_path.mkdir(parents=True, exist_ok=True)

    logger.info("Generating plots...")
    plot_dict = plots(found_models, out_path / PLOTS_DIR)

    logger.info("Comparing output models...")
    programs_run = []
    cif_models = [
        cif_file
        for cif_list in get_all_cif_files(found_models).values()
        for cif_file in cif_list
    ]
    indicies = get_gap_indicies(*cif_models)
    index_counter = 0

    alphafold_models = {"models": []}
    # TODO: Add AlphaFold 3 post-processing when needed

    boltz_models = {"models": []}
    if boltz_results_dir is not None:
        programs_run.append("Boltz")
        logger.info("Post-processing Boltz models...")
        for seed in bo.output.keys():
            for idx in bo.output[seed].keys():
                model = bo.output[seed][idx]["cif"]
                model.check_clashes()
                score_file = bo.output[seed][idx]["json"]
                plddt = model.residue_plddts
                if len(indicies) > 0:
                    plddt = insert_none_by_minus_one(indicies[index_counter], plddt)
                index_counter += 1
                model_data = get_model_data(
                    model,
                    plot_dict,
                    "Boltz",
                    plddt,
                    score_file,
                    out_path,
                )
                boltz_models["models"].append(model_data)

    chai_models = {"models": []}
    if chai_results_dir is not None:
        programs_run.append("Chai-1")
        logger.info("Post-processing Chai models...")
        for seed in co.output.keys():
            for idx in co.output[seed].keys():
                if idx >= 0:
                    model = co.output[seed][idx]["cif"]
                    model.check_clashes()
                    score_file = co.output[seed][idx]["scores"]
                    plddt = model.residue_plddts
                    if len(indicies) > 0:
                        plddt = insert_none_by_minus_one(indicies[index_counter], plddt)
                    index_counter += 1
                    model_data = get_model_data(
                        model,
                        plot_dict,
                        "Chai-1",
                        plddt,
                        score_file,
                        out_path,
                    )
                    chai_models["models"].append(model_data)

    combined_models = (
        alphafold_models["models"] + boltz_models["models"] + chai_models["models"]
    )
    # Generate output page
    import json
    import shutil

    from abcfold.abcfold import HTML_DIR, HTML_TEMPLATE
    from abcfold.html.html_utils import get_model_sequence_data, render_template
    from abcfold.output.file_handlers import superpose_models as superimpose_models

    (out_path / "output_models").mkdir(exist_ok=True)
    output_models = []
    logger.info("Preparing output model files...")
    for model in combined_models:
        cif_file = out_path.joinpath(model["model_path"])
        # if model["model_source"] == "AlphaFold3":
        #     output_name = "af3_model_" + model["model_id"][-1] + ".cif"
        # elif model["model_source"] == "Boltz":
        #     output_name = "boltz_model_" + model["model_id"][-1] + ".cif"
        # elif model["model_source"] == "Chai-1":
        #     output_name = "chai_model_" + model["model_id"][-1] + ".cif"
        output_name = f"{model['model_id']}.cif"
        output_model_path = out_path.joinpath("output_models").joinpath(output_name)
        shutil.copyfile(cif_file, output_model_path)
        output_models.append(output_model_path)
    if len(output_models) > 1:
        superimpose_models(output_models)

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
        "plotly_path": Path(plot_dict["plddt"]).relative_to(out_path).as_posix(),
        "chain_data": chain_data,
    }
    results_json = json.dumps(results_dict)

    if not out_path.joinpath(".feature_viewer").exists():
        shutil.copytree(HTML_DIR, out_path / ".feature_viewer")

    if len(programs_run) > 1:
        programs = (
            "Structure predictions for: "
            + ", ".join(programs_run[:-1])
            + " and "
            + programs_run[-1]
        )
    else:
        programs = "Structure predictions for: " + programs_run[0]

    # Create the index page
    logger.info("Generating output HTML page...")
    HTML_OUT = out_path.joinpath("index.html")
    html_out = Path(HTML_OUT).resolve()
    render_template(
        HTML_TEMPLATE,
        html_out,
        # kwargs appear as variables in the template
        abcfold_html_dir=".feature_viewer",
        programs=programs,
        results_json=results_json,
        version=0.1,
    )
    logger.info(f"Output page written to {HTML_OUT}")


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

        # The directory pattern is boltz_results_{run_id}_seed-{seed}
        self.name = (
            self.output_dirs[0].name.split("boltz_results_")[1].split("_seed-")[0]
        )

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
                pae = Af3Pae.from_boltz(pae_file.data, cif_file)
                pae_path = pae_file.pathway
                out_name = pae_path.with_name(f"{pae_path.stem}-af3.json")
                pae.to_file(out_name)

                if seed not in new_pae_files:
                    new_pae_files[seed] = []
                new_pae_files[seed].append(ConfidenceJsonFile(out_name))

        self.output = {
            seed: {
                i: {
                    "cif": cif_file,
                    "af3_pae": new_pae_files[seed][i],
                    "json": self.output[seed][i]["json"],
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

        # The directory pattern is chai_{run_id}_seed-{seed}
        self.name = self.output_dirs[0].name.split("chai_")[1].split("_seed-")[0]

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
        from abcfold.output.file_handlers import CifFile, FileTypes, NpyFile, NpzFile

        file_groups = {}

        for pathway in tqdm(self.output_dirs, desc="Processing Chai output"):
            seed = pathway.name.split("_")[-1]
            if seed not in file_groups:
                file_groups[seed] = {}

            # Handle PAE scores
            if (pae_file := pathway / "pae_scores.npy").exists():
                if -1 not in file_groups[seed]:
                    file_groups[seed][-1] = [NpyFile(str(pae_file))]
                else:
                    file_groups[seed][-1].append(NpyFile(str(pae_file)))

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
                        file_.pathway = (
                            file_.pathway.parent / f"stripped-{file_.pathway.name}"
                        )
                        file_.to_file(file_.pathway)
                        intermediate_dict["cif"] = file_

                if pae_file_data is not None:
                    # new_pae_path = (
                    #     file_.pathway.parent / f"pae_scores_model_{model_number}.npy"
                    # )
                    # if not new_pae_path.exists():
                    #     import shutil

                    #     shutil.copy(pae_file_data.pathway, new_pae_path)
                    # intermediate_dict["pae"] = NpyFile(str(new_pae_path))
                    intermediate_dict["pae"] = pae_file_data

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
                pae = Af3Pae.from_chai1(pae_file.data[i], cif_file)
                pae_path = pae_file.pathway
                out_name = pae_path.with_name(f"{pae_path.stem}-af3.json")
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
