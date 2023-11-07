import dash_bio
import atexit
import flask
import os
import threading

from dash import dcc, html
from dash.dependencies import Input, Output, State
from dash_app_func import DASH_CM_APP, APP_LAYOUT

core_path = os.getcwd()
if not os.path.exists(os.path.join(core_path, "tmp")):
    os.makedirs(os.path.join(core_path, "tmp"))

text_style = {
    "color": "#2E3138",
    "font-family": "Helvetica Neue",
    "width": "80%",
    "margin": "0 auto",
}

dash_cm_app = DASH_CM_APP(core_path=core_path)
# if dash_cm_app.get_temp_path() is None:
#     temp_path = tempfile.mkdtemp(
#     dir=os.path.join(self.core_path, "tmp")
# )
# atexit.register(dash_cm_app.clear_temp_directory)  # , dash_cm_app.get_temp_path()

app_layout = APP_LAYOUT(core_path=core_path)

igv_panel = dash_cm_app.core_igv_tab()
igv_panel.extend(
    [
        html.Div(id="selected-options"),
        dcc.Loading(parent_className="dashbio-loading", id="igv-output"),
        dcc.Store(id="igv-tracks-store", data={"tracks": []}),
    ]
)


def callbacks(_app):
    @_app.callback(
        Output("igv-output", "children"),
        [
            Input("method-dropdown", "value"),
            Input("cell_type-dropdown", "value"),
            Input("cm_qtls-switch", "checked"),
            Input("upload-data", "contents"),
            State("upload-data", "filename"),
            State("igv-tracks-store", "data"),
        ],
    )
    def update_selected_options(
        methods, cell_types, add_cmQTLs, file_content, file_name, tracks_state
    ):
        if (not methods or not cell_types) and (not file_content or not file_name):
            return "No options selected."
        else:
            updated_tracks = dash_cm_app.get_updated_igv_tracks(
                methods,
                cell_types,
                add_cmQTLs=add_cmQTLs,
                file_content=file_content,
                file_name=file_name,
                tracks_state=tracks_state,
            )
            return html.Div(
                children=[
                    dash_bio.Igv(
                        id="igv-chart",
                        genome="hg19",
                        reference=None,
                        tracks=updated_tracks,
                        # locus=["chr1:25500129-25800820"],
                    )
                ]
            )
            # return dash_cm_app.update_igv_tracks(
            #     methods,
            #     cell_types,
            #     add_cmQTLs=add_cmQTLs,
            #     file_content=file_content,
            #     file_name=file_name,
            #     tracks_state=tracks_state,
            # )

    # @_app.callback(Output("igv-chart", "tracks"), Input("igv-tracks", "data"))
    # def update_igv_chart_tracks(updated_tracks):
    #     return updated_tracks


app = app_layout.run_standalone_app(dash_cm_app.layout(igv_panel), callbacks)
server = app.server

if core_path.endswith("/"):
    core_path = core_path[:-1]


@server.route(core_path + "/<path:path>")
def send_bw(path):
    return flask.send_from_directory(core_path, path)


if __name__ == "__main__":
    app.run_server(debug=True, port=8050, host="0.0.0.0")  # mode="jupyterlab",
