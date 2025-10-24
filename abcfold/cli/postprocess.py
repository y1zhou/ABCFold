"""Run the post-processing steps for ABCFold."""

import logging
from pathlib import Path
from typing import Annotated

import typer

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

        boltz_out_dirs = list(boltz_results_dir.glob("boltz_seed_*/boltz_results*"))
        run_id = boltz_out_dirs[0].name.split("boltz_results_")[1]
        bo = PatchedBoltzOutput(boltz_out_dirs, None, run_id, False)
        found_models.append(bo)

    if chai_results_dir is not None:
        logger.info("Post-processing Chai results...")

        chai_out_dirs = list(chai_results_dir.glob("chai_seed_*"))
        run_seed = chai_out_dirs[0].name
        run_log_path = next((chai_out_dirs[0].parent / "logs").glob("*.log"))
        run_id = run_log_path.name.split(run_seed)[0]
        co = PatchedChaiOutput(chai_out_dirs, None, run_id, False)
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
    for model in combined_models:
        cif_file = out_path.joinpath(model["model_path"])
        if model["model_source"] == "AlphaFold3":
            output_name = "af3_model_" + model["model_id"][-1] + ".cif"
        elif model["model_source"] == "Boltz":
            output_name = "boltz_model_" + model["model_id"][-1] + ".cif"
        elif model["model_source"] == "Chai-1":
            output_name = "chai_model_" + model["model_id"][-1] + ".cif"
        shutil.copy(
            cif_file,
            out_path.joinpath("output_models").joinpath(output_name),
        )
        output_models.append(out_path.joinpath("output_models").joinpath(output_name))
    if len(output_models) > 1:
        superimpose_models(output_models)

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
):
    """Serve the output HTML page with a local web server."""
    import os
    import socketserver
    import sys
    import webbrowser

    from abcfold.abcfold import PORT
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
    output_open_html_script("open_output.py", port=PORT)

    try:
        # Start the server
        with socketserver.TCPServer(("", PORT), NoCacheHTTPRequestHandler) as httpd:
            logger.info(f"Serving at port {PORT}: http://localhost:{PORT}/index.html")
            logger.info("Press Ctrl+C to stop the server")
            # Open the main HTML page in the default web browser
            webbrowser.open(f"http://localhost:{PORT}/index.html")
            # Keep the server running
            httpd.serve_forever()
    except KeyboardInterrupt:
        logger.info("Server stopped")
        httpd.server_close()
        sys.exit(0)


class PatchedBoltzOutput(BoltzOutput):
    """Patched BoltzOutput to bypass get_input_fasta()."""

    def __init__(
        self,
        boltz_output_dirs: list[str | Path],
        input_params: dict,
        name: str,
        save_input: bool = False,
    ):
        """Process the output of a Boltz run.

        Ref: `abcfold.abcfold.output.boltz:BoltzOutput`
        """
        self.output_dirs = [Path(x) for x in boltz_output_dirs]
        self.input_params = input_params
        self.name = name
        self.save_input = save_input

        parent_dir = self.output_dirs[0].parent
        new_parent = parent_dir / f"boltz_{self.name}"
        new_parent.mkdir(parents=True, exist_ok=True)

        if self.save_input:
            boltz_yaml = list(parent_dir.glob("*.yaml"))[0]
            if boltz_yaml.exists():
                boltz_yaml.rename(new_parent / "boltz_input.yaml")

            boltz_msas = list(parent_dir.glob("*.a3m"))
            if boltz_msas:
                for boltz_msa in boltz_msas:
                    if boltz_msa.exists():
                        boltz_msa.rename(new_parent / boltz_msa.name)

        new_output_dirs = []
        for output_dir in self.output_dirs:
            if output_dir.name.startswith("boltz_results_"):
                new_path = new_parent / output_dir.name
                output_dir.rename(new_path)
                new_output_dirs.append(new_path)
            else:
                new_output_dirs.append(output_dir)
        self.output_dirs = new_output_dirs

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


class PatchedChaiOutput(ChaiOutput):
    """Patched ChaiOutput to bypass get_input_fasta()."""

    def get_input_fasta(self):
        """Return None to bypass input fasta retrieval."""
        return
