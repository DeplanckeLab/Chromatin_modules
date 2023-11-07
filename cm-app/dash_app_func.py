import base64
import dash
import dash_bio
import dash_mantine_components as dmc
import os
import tempfile
import shutil
from dash import dcc, html

# from jupyter_dash import JupyterDash


class DASH_CM_APP:
    global_temp_path = None

    def __init__(self, core_path):
        self.core_path = core_path

    @staticmethod
    def core_igv_tab():
        CM_mapping_methods = [
            {"value": "VCM", "label": "VCMtools"},
            {"value": "CRD", "label": "Clomics"},
            {"value": "PHM", "label": "PHM"},
        ]
        cell_types = [
            {"value": "LCL", "label": "LCL"},
            {"value": "FIB", "label": "FIB"},
            {"value": "Monocytes", "label": "Monocytes"},
            {"value": "Neutrophils", "label": "Neutrophils"},
            {"value": "Tcells", "label": "Tcells"},
        ]
        return [
            html.Div(
                id="genome",
                className="fullwidth-app-controls-name",
                children=html.Div(
                    [
                        html.A(
                            "Genome: ",
                            style={
                                "font-weight": "bold",
                            },
                        ),
                        "Human (GRCh37/hg19)",
                    ]
                ),
            ),
            html.Div(
                id="method",
                className="fullwidth-app-controls-name",
                children="Select CM mapping strategy:",
                style={
                    "font-weight": "bold",
                },
            ),
            dcc.Dropdown(
                id="method-dropdown",
                options=[
                    {
                        "label": method_dict.get("label"),
                        "value": method_dict.get("value"),
                    }
                    for method_dict in CM_mapping_methods
                ],
                multi=True,
            ),
            html.Div(
                id="cell_type",
                className="fullwidth-app-controls-name",
                children="Select cell type:",
                style={
                    "font-weight": "bold",
                },
            ),
            dcc.Dropdown(
                id="cell_type-dropdown",
                options=[
                    {
                        "label": cell_type_dict.get("label"),
                        "value": cell_type_dict.get("value"),
                    }
                    for cell_type_dict in cell_types
                ],
                multi=True,
            ),
            html.Div(
                id="cm_qtls",
                className="fullwidth-app-controls-name",
                children="Add cmQTLs:",
                style={
                    "font-weight": "bold",
                },
            ),
            dmc.Switch(
                id="cm_qtls-switch",
                size="sm",
                checked=False,
                color="#45649D",
            ),
            dcc.Upload(
                id="upload-data",
                children=html.Div(
                    [
                        "Drag and Drop or ",
                        html.A(
                            "Select Files",
                            style={"font-weight": "bold", "color": "#45649D"},
                        ),
                    ]
                ),
                style={
                    "width": "100%",
                    "height": "60px",
                    "lineHeight": "60px",
                    "borderWidth": "1px",
                    "borderStyle": "dashed",
                    "borderRadius": "5px",
                    "textAlign": "center",
                    "margin": "10px",
                },
                multiple=True,
                max_size=500 * 1024 * 1024,  # 500 MB (adjust as needed)
            ),
        ]

    # @staticmethod
    # def header_colors():
    #     return {
    #         "bg_color": "#2E3138",
    #         "font_color": "white",
    #     }

    def layout(self, igv_panel):
        cm_intro_text = """
            _**Chromatin modules**_ (CMs) represent sub-TAD hubs encompassing interactions
            between cis-regulatory elements, such as promoters and enhancers.
            CMs are charachterized by
            - interaction of active promoters with one or more active enhancers simultaneously
            (median of two-four enhancers),
            - preferential scaling at a sub-TAD level and therefore often being <100 kb
            in size (with median sizes ranging from 32 to 70 kb),
            - physical proximity of interacting elements within modules that show
            enrichment in 3D contacts,
            - deposition of histone modifications and binding of TFs in such regions
            in a coordinated fashion with transcriptional
            activity being a strong predictor of module formation.

            For more details, please refer to [van Mierlo, Pushkarev,
            Trends in Genetics, 2023](https://www.cell.com/trends/genetics/fulltext/S0168-9525(22)00290-6)).
        """
        cm_mapping_intro = """
            In [Pushkarev, van Mierlo et al., bioRxiv, 2023](https://www.biorxiv.org/content/10.1101/2023.10.11.561870v1)
            we evaluated three computational methods for CM mapping based on epigenome data
            obtained for a range of individuals.
        """
        cm_mapping_text_pt1 = """
            In this web portal, we provide the CM browser that displays chromatin modules
            mapped using H3K27ac and H3K4me1 ChIP-seq data from a large number of individuals
            for five different cell types:
        """
        cm_mapping_text_pt2 = """
            and three computational approaches:
            - [VCMtools](https://doi.org/10.1016/j.cell.2015.08.001),
            - [Clomics](https://www.science.org/doi/10.1126/science.aat8266),
            - [PHM](https://www.nature.com/articles/s41588-018-0278-6).

            CM are visualized as stretches of linked rectangles, where each rectangle represents
            a ChIP-seq peak for either H3K27ac or H3K4me1. CMs are colored with respect to the cell type of choice.
            Each CM is described by its chromosome positioning, start (the leftmost CM peak)
            and end (the rightmost CM peak) coordinates. The "score" field indicates average
            correlation strength of CM peaks.

            To visualize CMs:

            1. Select a CM mapping method (**Select CM mapping strategy**),
            2. Select a cell type of interest (**Select cell type**).

            We also provide an option to visualize the genetic variants with the smallest
            p-value association to the mapped CM (**cmQTLs**).
            Note that only significant variants are indicated, and in cases where multiple variants had the
            same p-value of association with CM activity, only one variant is shown.
        """
        page_1 = dcc.Tab(
            label="About Chromatin Modules",
            value="definition",
            children=html.Div(
                className="control-tab",
                children=[
                    html.H3(
                        className="definition",
                        children="Chromatin modules and their definition",
                    ),
                    dcc.Markdown(cm_intro_text),
                    html.Img(
                        src=os.path.join(self.core_path, "assets", "Figure2_v1.png"),
                        alt="Chromatin modules",
                        style={
                            "max-width": "90%",
                            "height": "auto",
                            "display": "block",
                            "margin": "0 auto",
                            "margin-bottom": "10px",
                        },
                    ),
                    html.Hr(),
                    html.H3(
                        className="cm_mapping",
                        children="Mapping of CMs across cell types",
                    ),
                    dcc.Markdown(cm_mapping_intro),
                    html.Img(
                        src=os.path.join(
                            self.core_path,
                            "assets",
                            "Figure1_v1.png",
                        ),
                        alt="Chromatin modules",
                        style={
                            "max-width": "90%",
                            "height": "auto",
                            "display": "block",
                            "margin": "0 auto",
                            "margin-bottom": "10px",
                        },
                    ),
                    html.Hr(),
                    html.H3(
                        className="functionality",
                        children="Functionality of the browser",
                    ),
                    dcc.Markdown(cm_mapping_text_pt1),
                    html.Img(
                        src=os.path.join(
                            self.core_path,
                            "assets",
                            "cell_types.png",
                        ),
                        alt="cell types",
                        style={
                            "float": "left",
                            "width": "11%",
                            "height": "auto",
                            "display": "block",
                            "margin-right": "1500px",
                            "margin-bottom": "10px",
                        },
                    ),
                    dcc.Markdown(cm_mapping_text_pt2),
                    html.Hr(),
                ],
            ),
        )
        page_2 = dcc.Tab(
            label="Explore Chromatin Modules",
            value="explore-cms",
            children=html.Div(
                className="control-tab",
                children=igv_panel,
            ),
        )
        page_3 = dcc.Tab(
            label="Download data",
            value="download_cms",
            children=html.Div(
                className="control-tab", children=[dcc.Markdown("In progress...")]
            ),
        )
        return html.Div(
            id="igv-body",
            className="app-body",
            children=[
                html.Div(
                    style={
                        "color": "#2E3138",
                        "font-family": "Helvetica Neue",
                        "width": "80%",
                        "margin": "0 auto",
                    },
                    id="igv-control-tabs",
                    className="control-tabs",
                    children=[
                        dcc.Tabs(
                            id="igv-tabs",
                            value="definition",
                            children=[page_1, page_2, page_3],
                        )
                    ],
                ),
            ],
        )

    def get_updated_igv_tracks(
        self,
        methods=None,
        cell_types=None,
        add_cmQTLs=None,
        file_content=None,
        file_name=None,
        tracks_state=None,
    ):
        updated_tracks = tracks_state.get("tracks", [])
        if (methods is not None) and (cell_types is not None):
            for method in methods:
                for cell_type in cell_types:
                    if method == "VCM":
                        method_name = "VCMtools"
                        suffix = "_0.001"
                    elif method == "CRD":
                        method_name = "Clomics"
                        suffix = ""
                    else:
                        method_name = "PHM"
                        suffix = "_0.8"
                    updated_tracks.append(
                        {
                            "name": method_name + ", " + cell_type,
                            "url": os.path.join(
                                self.core_path,
                                "data",
                                "CMs_colored_by_cell_type",
                                method_name,
                                "_".join(
                                    [
                                        cell_type,
                                        method + suffix,
                                        "colored_by_cell_type.tracks.bed",
                                    ]
                                ),
                            ),
                            "displayMode": "EXPANDED",
                            "type": "annotation",
                            "format": "bed",
                            "removable": True,
                        }
                    )
                    if add_cmQTLs:
                        updated_tracks.append(
                            {
                                "name": method_name + ", " + cell_type + ", cmQTLs",
                                "url": os.path.join(
                                    self.core_path,
                                    "data",
                                    "cmQTLs",
                                    method_name,
                                    "_".join([cell_type, method, "cmQTLs.bed"]),
                                ),
                                "displayMode": "EXPANDED",
                                "type": "annotation",
                                "format": "bed",
                                "removable": True,
                            }
                        )
        if file_content is not None:
            for fc, fn in zip(file_content, file_name):
                DASH_CM_APP.global_temp_path = tempfile.mkdtemp(
                    dir=os.path.join(self.core_path, "tmp")
                )
                _, content_string = fc.split(",")
                decoded_data = base64.b64decode(content_string)
                temp_path = os.path.join(DASH_CM_APP.global_temp_path, fn)
                with open(temp_path, "wb") as temp_file:
                    temp_file.write(decoded_data)
                updated_tracks.append(
                    {
                        "name": ".".join(fn.split(".")[:-1]),
                        "url": temp_path,
                        "removable": True,
                    }
                )
        return updated_tracks
        # return html.Div(
        #     children=[
        #         dash_bio.Igv(
        #             id="igv-chart",
        #             genome="hg19",
        #             reference=None,
        #             tracks=tracks,
        #             # locus=["chr1:25500129-25800820"],
        #         )
        #     ]
        # )

    # def get_temp_path(self):
    #     return DASH_CM_APP.global_temp_path

    # def clear_temp_directory(self):
    #     if (
    #         os.path.exists(DASH_CM_APP.global_temp_path)
    #         and "DASH" in DASH_CM_APP.global_temp_path
    #     ):
    #         shutil.rmtree(DASH_CM_APP.global_temp_path)


class APP_LAYOUT:
    def __init__(self, core_path):
        self.core_path = core_path

    def run_standalone_app(self, layout, callbacks):
        app = dash.Dash(
            __name__,
            external_stylesheets=[
                os.path.join(self.core_path, "footer_style.css"),
            ],
        )
        # app = JupyterDash(
        #     __name__,
        #     external_stylesheets=[
        #         os.path.join(self.core_path, "footer_style.css"),
        #     ],
        # )

        # Handle callback to component with id "fullband-switch"
        app.config["suppress_callback_exceptions"] = True

        # Assign layout
        app.layout = self.app_page_layout(
            page_layout=layout,
            light_logo=True,
            bg_color="dimgray",
            font_color="#F3F6FA",
        )

        # Register all callbacks
        callbacks(app)

        # return app object
        return app

    def app_page_layout(self, page_layout, light_logo, bg_color, font_color):
        return html.Div(
            id="main_page",
            children=[
                dcc.Location(id="url", refresh=False),
                html.Div(
                    id="app-page-header",
                    children=[
                        html.H1(
                            "Chromatin Modules",
                            style={
                                "display": "inline",
                                "font-family": "Helvetica Neue",
                                "color": "white",
                                "margin-left": "10px",
                            },
                        ),
                        html.A(
                            html.Img(
                                src="data:image/png;base64,{}".format(
                                    base64.b64encode(
                                        open(
                                            os.path.join(
                                                self.core_path,
                                                "assets",
                                                "github-mark-white.png",
                                            ).format("Light-" if light_logo else ""),
                                            "rb",
                                        ).read()
                                    ).decode()
                                ),
                                style={
                                    "width": "30px",
                                    "height": "30px",
                                    "position": "absolute",
                                    "top": "5px",
                                    "right": "10px",
                                },
                            ),
                            href="https://github.com/DeplanckeLab/Chromatin_modules",
                        ),
                    ],
                    style={
                        "background": bg_color,
                        "font-family": "Helvetica Neue",
                        "color": font_color,
                    },
                ),
                html.Div(id="app-page-content", children=page_layout),
                html.Div(
                    [
                        html.P(
                            [
                                "Please contact ",
                                html.A(
                                    "Olga Pushkarev",
                                    href="mailto:olga.pushkareva@epfl.ch",
                                ),
                                " or ",
                                html.A(
                                    "Guido van Mierlo",
                                    href="mailto:guido.vanmierlo@epfl.ch",
                                ),
                                " if you have any questions or suggestions regarding the portal.",
                            ],
                            style={
                                "font-family": "Helvetica Neue",
                                "color": "#474747",
                            },
                        ),
                    ],
                    className="footer",
                ),
            ],
        )
