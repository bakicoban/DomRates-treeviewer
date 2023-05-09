# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc

from dash.dependencies import Input, Output
from utils import (
    create_paths_file,
    create_tree,
)

#app = dash.Dash(__name__)
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

app.title = "Phylogeny Tree Explorer"
server = app.server

protein_name = "test"
proteins = ["test", "test2"]
tree_fig = {}

tree_file, metadata_file = create_paths_file(protein_name)

fig = create_tree(protein_name, tree_file, metadata_file)
tree_fig[tree_file] = fig


######################################### MAIN APP #########################################
app.layout = html.Div(
    [
        # Banner display
        html.Div(
            className="header-title",
            children=[
                html.H2(
                    id="title",
                    children="Phylogeny trees of DomRates results",
                ),
                html.Div(
                    id="learn_more",
                    children=[
                        html.Img(className="logo", src=app.get_asset_url("logo.png"))
                    ],
                ),
            ],
        ),
        html.Div(
            id="row",
            children=[
                html.Div(
                    id="controls",
                    className="row div-row div-card",
                    children=[
                        html.Div(
                            id="dataset-picker",
                            children=[
                                html.Div(
                                    className="six columns",
                                    children=[
                                        html.H6(children="Protein name"),
                                        dcc.Dropdown(
                                            id="d_protein-name",
                                            options=[
                                                {
                                                    "label": proteins[i],
                                                    "value": proteins[i],
                                                }
                                                for i in range(len(proteins))
                                            ],
                                            value="test",
                                        ),
                                        html.Div(id="output-container"),
                                    ],
                                ),
                            ],
                        ),

                    ],
                ),
                dcc.Graph(id="phylogeny-graph", className="div-card", figure=fig),
            ],
        ),
    ]
)


"""app.layout = html.Div(
    [
        # Banner display
        html.Div(
            className="header-title",
            children=[
                html.H2(
                    id="title",
                    children="Phylogeny trees of DomRates results",
                ),
                html.Div(
                    id="learn_more",
                    children=[
                        html.Img(className="logo", src=app.get_asset_url("logo.png"))
                    ],
                ),
            ],
        ),
        html.Div(
            id="row",
            children=[
                html.Div(
                    id="controls",
                    className="row div-row div-card",
                    children=[
                        html.Div(
                            id="dataset-picker",
                            children=[
                                html.Div(
                                    className="six columns",
                                    children=[
                                        html.H6(children="Protein name"),
                                        dcc.Dropdown(
                                            id="d_protein-name",
                                            options=[],
                                            placeholder="Type something...",
                                            value=None,
                                        ),
                                        html.Div(id="output-container"),
                                    ],
                                ),
                            ],
                        ),
                    ],
                ),
                dcc.Graph(id="phylogeny-graph", className="div-card", figure=fig),
            ],
        ),
    ]
)"""




######################################### UPDATING FIGURES #########################################
@app.callback(Output("output-container", "children"), [Input("d_protein-name", "value")])
def _update_legend_gene(protein_name):
    return 'You have selected "{}" protein'.format(protein_name)



@app.callback(
    Output("phylogeny-graph", "figure"),
    [
        Input("d_protein-name", "value"),
    ],
)



def update_phylogeny_tree(protein_name):
    protein_name = protein_name.lower()
    tree_file_filtred, metadata_file_filtred = create_paths_file(protein_name)

    if tree_file_filtred in tree_fig:
        fig = tree_fig[tree_file_filtred]
    else:
        fig = create_tree(protein_name, tree_file_filtred, metadata_file_filtred)
        tree_fig[tree_file_filtred] = fig
    return fig

# Running the server
if __name__ == "__main__":
    app.run_server(debug=True)
