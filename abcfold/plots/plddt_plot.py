"""Plotting functions for pLDDT distributions."""

import logging
from pathlib import Path
from typing import Any

import plotly.express as px
import plotly.graph_objects as go
import plotly.offline as pyo

from abcfold.output.file_handlers import CifFile
from abcfold.output.utils import get_gap_indicies, insert_none_by_minus_one

logger = logging.getLogger("logger")


def plot_plddt(
    cif_models_dict: dict[str, list[CifFile]],
    output_name: str | Path,
    line_width: float = 1.6,
    dash: str = "dot",
    show: bool = False,
    chain_line_occupancy: float = 0.9,
    include_plotlyjs: bool = True,
) -> None:
    """Plot the pLDDT distribution of the models in the dictionary of cif models.

    Args:
        cif_models_dict: Dictionary of cif models to plot. The keys are the source of
            the models and the values are lists of CifFile objects.
            e.g. {"Alphafold3": [CifFile, CifFile, ...], "Boltz": [CifFile, ...],
            "Chai-1": [CifFile, ...]}

        output_name: Path to the output html file.
        line_width: Width of the lines in the plot.
        dash: Dash style of the lines in the plot.
        show: If True, the plot will be displayed in the browser.
        chain_line_occupancy: Opacity of the vertical lines that separate the chains.
        include_plotlyjs: If True, the plotly.js library will be included in the html file.

    Returns:
        None

    Outputs:
        An html file with the plot.

    """
    fig = go.Figure()
    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor="LightGrey")
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor="LightGrey")

    colours = list(px.colors.qualitative.T10)
    method_colours = {
        "Alphafold3": px.colors.qualitative.Set1,
        "Boltz": px.colors.qualitative.Set2,
        "Chai-1": px.colors.qualitative.Prism,
    }

    line_ranges: dict = {}

    cif_models = [
        cif_file for cif_files in cif_models_dict.values() for cif_file in cif_files
    ]
    indicies = get_gap_indicies(*cif_models)
    indicies_index = 0

    added_model_order: list[tuple[str, int, int, str, go.Scatter]] = []
    for method, cif_models in cif_models_dict.items():
        for cif_model in sorted(cif_models, key=lambda f: f.name):
            cif_model_name_parts = cif_model.name.split("_")
            model_index = int(cif_model_name_parts[-1])
            model_seed = int(cif_model_name_parts[-2].split("-")[-1])
            model_name = f"Seed {model_seed}-M{model_index}"

            color_list = method_colours.get(method, colours)
            color = color_list[model_index % len(color_list)]

            plddt = cif_model.residue_plddts

            if len(indicies) > 0:
                plddt = insert_none_by_minus_one(indicies[indicies_index], plddt)

            indicies_index += 1
            chain_ranges = {
                chain: len(plddt)
                for chain, plddt in cif_model.get_plddt_per_residue().items()
            }
            line_ranges = {
                chain: max(chain_ranges[chain], line_ranges.get(chain, 0))
                for chain in chain_ranges
            }

            trace = go.Scatter(
                x=list(range(len(plddt))),
                y=plddt,
                mode="lines",
                legendgroup=method,
                legendgrouptitle_text=Bold(method),
                name=model_name,
                line=dict(dash=dash, width=line_width, color=color),
                visible=True,  # Ensure traces start as visible
                showlegend=True,
            )
            added_model_order.append(
                (method, model_seed, model_index, model_name, trace)
            )

    # Order lines by method, seed, and model index
    added_model_order = sorted(added_model_order, key=lambda x: (x[0], x[1], x[2]))
    for *_, trace in added_model_order:
        fig.add_trace(trace)

    fig.update_layout(
        legend=dict(
            title=dict(
                text="Click to Show/Hide Methods",
                font=dict(size=12, color="black"),
            )
        )
    )

    counter = 0
    colour_index = 0
    for chain, chain_range in line_ranges.items():
        counter += chain_range
        chain_name = f"Chain {chain}"

        fig.add_vline(
            x=counter - 1,
            line=dict(color=colours[colour_index % len(colours)], dash="dash"),
            opacity=chain_line_occupancy,
            annotation_text=Bold(chain_name),
            annotation_font_size=15,
            annotation_position="top left",
            annotation_textangle=-90,
        )

        colour_index += 1

    # Create buttons for each model
    buttons = []

    # Add a button to show all traces
    buttons.append(
        dict(
            method="update",
            args=[{"visible": [True] * len(fig.data)}, {"showlegend": True}],
            label="All",
        )
    )

    # Add buttons for each individual model
    added_model_names = []
    for _, _, _, model_name, _ in added_model_order:
        if model_name not in added_model_names:
            added_model_names.append(model_name)

    for model_name in added_model_names:
        button: dict[str, Any] = dict(
            method="update",
            args=[
                {"visible": [model_name == n for _, _, _, n, _ in added_model_order]},
                {"showlegend": True},
            ],
            label=model_name,
        )
        buttons.append(button)

    # Add the updatemenu to the layout
    fig.update_layout(
        updatemenus=[dict(showactive=True, buttons=buttons)],
        xaxis_title=Bold("Residue Number"),
        yaxis_title=Bold("pLDDT Score"),
        title=Bold("pLDDT Distribution"),
        plot_bgcolor="white",
    )

    if show:
        fig.show()

    output_name = Path(output_name)

    if output_name.suffix == "":
        output_name = output_name.with_suffix(".html")

    if include_plotlyjs:
        fig.write_html(str(output_name))
    else:
        div = pyo.plot(fig, include_plotlyjs=False, output_type="div")
        output_name = output_name.with_suffix(".div.html")
        with open(output_name, "w") as f:
            f.write(div)


def Bold(string):
    """Make a string bold in plotly annotations."""
    return f"<b>{string}</b>"
