import dash
from dash import Input, Output, State, dcc, html
import plotly.express as px
import numpy as np
import pyvista as pv
from dash_vtk.utils import to_mesh_state
import dash_vtk
import os

# Import local libs
from lib.data_loader import load_sim_snapshot, get_volumetric_averages, get_available_timesteps, discover_projects

# 1. Update Project List
@dash.callback(
    [Output('pressure-run-selector', 'options'),
     Output('pressure-run-selector', 'value')],
    Input('pressure-btn-refresh', 'n_clicks')
)
def update_project_list(n):
    prefixes = discover_projects("exports")
    options = [{'label': p, 'value': p} for p in prefixes]
    default_val = options[0]['value'] if options else None
    return options, default_val

# 2. Update Time Slider
@dash.callback(
    [Output('pressure-time-slider', 'max'),
     Output('pressure-time-slider', 'marks'),
     Output('pressure-time-slider', 'value')],
    Input('pressure-run-selector', 'value')
)
def update_time_slider(prefix):
    steps = get_available_timesteps(prefix)
    if not steps:
        return 0, {0: '0'}, 0
    marks = {s: str(s) for s in steps[::max(1, len(steps)//10)]}
    return max(steps), marks, min(steps)

# 3. Main Rendering Logic
@dash.callback(
    [Output('pressure-tab-content', 'children'),
     Output('pressure-z-slider', 'max'),
     Output('pressure-z-slider', 'marks')],
    [Input('pressure-tabs', 'active_tab'),
     Input('pressure-run-selector', 'value'),
     Input('pressure-time-slider', 'value'),
     Input('pressure-var-selector', 'value'),
     Input('pressure-z-slider', 'value')]
)
def render_pressure_content(active_tab, prefix, t_step, var_name, z_layer):
    if not prefix:
        return html.Div("Please select a simulation project."), 0, {0: '0'}
    
    try:
        data = load_sim_snapshot(prefix, t_step)
        nz = data['grid']['nz']
        marks = {i: str(i) for i in range(nz)}
        
        # Color Scale: Viridis or Blue-Red for pressure
        color_scale = "Viridis"

        if active_tab == "tab-2d":
            field_key = var_name
            if field_key not in data['fields']:
                field_key = list(data['fields'].keys())[0]
            
            field = data['fields'][field_key]
            slice_data = field[min(z_layer, nz-1), :, :]
            
            fig = px.imshow(
                slice_data,
                color_continuous_scale=color_scale,
                labels=dict(x="X", y="Y", color=var_name),
                title=f"{var_name} Profile - Step {t_step} (Layer {z_layer})"
            )
            fig.update_layout(
                template="plotly_dark",
                paper_bgcolor="rgba(0,0,0,0)",
                plot_bgcolor="rgba(0,0,0,0)",
                margin=dict(l=20, r=20, t=40, b=20)
            )
            return dcc.Graph(figure=fig, style={"height": "600px"}), nz-1, marks

        elif active_tab == "tab-3d":
            grid = data['grid']
            spacing = data['spacing']
            
            vti_field_name = var_name
            if vti_field_name not in data['fields']:
                vti_field_name = list(data['fields'].keys())[0]
            
            field_data = data['fields'][vti_field_name].flatten(order='F')
            
            min_val, max_val = float(np.min(field_data)), float(np.max(field_data))
            if max_val == min_val: max_val += 1.0
            
            mesh = pv.ImageData(
                dimensions=(grid['nx'], grid['ny'], grid['nz']),
                spacing=(spacing['dx'], spacing['dy'], spacing['dz'])
            )
            mesh.point_data[vti_field_name] = field_data
            mesh.set_active_scalars(vti_field_name)
            surface = mesh.extract_surface()
            mesh_state = to_mesh_state(surface)
            
            vtk_view = dash_vtk.View(
                id="pressure-vtk-view",
                children=[
                    dash_vtk.GeometryRepresentation(
                        id="pressure-vtk-repr",
                        children=[dash_vtk.Mesh(id="pressure-vtk-mesh", state=mesh_state)],
                        colorMapPreset="Cool to Warm",
                        colorDataRange=[float(min_val), float(max_val)],
                        showCubeAxes=True
                    )
                ],
                cameraPosition=[1, 1, 1],
                background=[0.05, 0.05, 0.05]
            )
            return html.Div(vtk_view, id="pressure-vtk-container", style={"height": "100%", "width": "100%"}), nz-1, marks
            
        elif active_tab == "tab-trend":
            time_steps, trend_data = get_volumetric_averages(prefix)
            if not time_steps:
                return html.Div("No history discovered."), nz-1, marks
                
            y_data = trend_data.get(var_name)
            if y_data is None:
                y_data = list(trend_data.values())[0]

            fig = px.line(
                x=time_steps, y=y_data,
                labels={"x": "Time Step", "y": f"Avg {var_name}"},
                template="plotly_dark"
            )
            fig.add_vline(x=t_step, line_dash="dash", line_color="orange")
            fig.update_layout(
                paper_bgcolor="rgba(0,0,0,0)",
                plot_bgcolor="rgba(0,0,0,0)"
            )
            return dcc.Graph(figure=fig, style={"height": "600px"}), nz-1, marks
            
    except Exception as e:
        return html.Div(f"Extraction Error: {e}", className="text-danger p-4"), 0, {0: '0'}
