import dash
from dash import Input, Output, State, dcc, html
import plotly.express as px
import numpy as np
import pyvista as pv
from dash_vtk.utils import to_mesh_state
import dash_vtk
import glob
import os
import sys
from functools import lru_cache

# Import local libs
from lib.data_loader import load_sim_snapshot, get_volumetric_averages, get_available_timesteps, discover_projects

@lru_cache(maxsize=32)
def get_cached_trend(prefix):
    return get_volumetric_averages(prefix)

# 1. Update Project List
@dash.callback(
    [Output('run-selector', 'options'),
     Output('run-selector', 'value')],
    Input('btn-refresh', 'n_clicks')
)
def update_project_list(n):
    prefixes = discover_projects("exports")
    options = [{'label': p, 'value': p} for p in prefixes]
    default_val = options[0]['value'] if options else None
    return options, default_val

# 2. Update Time Slider
@dash.callback(
    [Output('time-slider', 'max'),
     Output('time-slider', 'marks'),
     Output('time-slider', 'value')],
    Input('run-selector', 'value')
)
def update_time_slider(prefix):
    steps = get_available_timesteps(prefix)
    if not steps:
        return 0, {0: '0'}, 0
        
    marks = {s: str(s) for s in steps[::max(1, len(steps)//10)]}
    return max(steps), marks, min(steps)

# 3. Main Rendering Logic
@dash.callback(
    [Output('tab-content', 'children'),
     Output('z-slider', 'max'),
     Output('z-slider', 'marks')],
    [Input('tabs', 'active_tab'),
     Input('run-selector', 'value'),
     Input('time-slider', 'value'),
     Input('var-selector', 'value'),
     Input('z-slider', 'value')]
)
def render_reservoir_content(active_tab, prefix, t_step, var_name, z_layer):
    if not prefix:
        return html.Div("Please select a simulation project."), 0, {0: '0'}
    
    try:
        data = load_sim_snapshot(prefix, t_step)
        nz = data['grid']['nz']
        marks = {i: str(i) for i in range(nz)}
        
        if active_tab == "tab-2d":
            field = data['fields'][var_name]
            slice_data = field[min(z_layer, nz-1), :, :]
            fig = px.imshow(
                slice_data,
                color_continuous_scale="Viridis",
                labels=dict(x="X", y="Y", color=var_name),
                title=f"{var_name} Profile - Step {t_step} (Layer {z_layer})"
            )
            fig.update_layout(template="plotly_dark")
            return dcc.Graph(figure=fig, style={"height": "600px"}), nz-1, marks

        elif active_tab == "tab-3d":
            grid = data['grid']
            spacing = data['spacing']
            field_data = data['fields'][var_name].flatten(order='F')
            
            min_val, max_val = float(np.min(field_data)), float(np.max(field_data))
            if max_val == min_val: max_val += 1.0
            
            mesh = pv.ImageData(
                dimensions=(grid['nx'], grid['ny'], grid['nz']),
                spacing=(spacing['dx'], spacing['dy'], spacing['dz'])
            )
            mesh.point_data[var_name] = field_data
            mesh.set_active_scalars(var_name)
            surface = mesh.extract_surface()
            mesh_state = to_mesh_state(surface)
            
            vtk_view = dash_vtk.View(
                id="vtk-view",
                children=[
                    dash_vtk.GeometryRepresentation(
                        id="vtk-repr",
                        children=[dash_vtk.Mesh(id="vtk-mesh", state=mesh_state)],
                        colorMapPreset="Cool to Warm",
                        colorDataRange=[float(min_val), float(max_val)],
                        showCubeAxes=True
                    )
                ],
                cameraPosition=[1, 1, 1],
                background=[0, 0, 0]
            )
            return html.Div(vtk_view, id="vtk-container", style={"height": "100%", "width": "100%"}), nz-1, marks
            
        elif active_tab == "tab-trend":
            time_steps, trend_data = get_cached_trend(prefix)
            if not time_steps:
                return html.Div("No history discovered."), nz-1, marks
                
            fig = px.line(
                x=time_steps, y=trend_data[var_name],
                labels={"x": "Time Step", "y": f"Avg {var_name}"},
                template="plotly_dark"
            )
            fig.add_vline(x=t_step, line_dash="dash", line_color="orange")
            return dcc.Graph(figure=fig, style={"height": "600px"}), nz-1, marks
            
    except Exception as e:
        return html.Div(f"Error: {e}"), 0, {0: '0'}
