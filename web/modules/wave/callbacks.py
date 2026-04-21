import dash
from dash import Input, Output, State, dcc, html
import plotly.express as px
import numpy as np
import os

# Import local libs
from lib.data_loader import load_sim_snapshot, get_volumetric_averages, get_available_timesteps, discover_projects

# 1. Update Project List
@dash.callback(
    [Output('wave-run-selector', 'options'),
     Output('wave-run-selector', 'value')],
    Input('wave-btn-refresh', 'n_clicks')
)
def update_project_list(n):
    prefixes = discover_projects("exports")
    options = [{'label': p, 'value': p} for p in prefixes]
    default_val = options[0]['value'] if options else None
    return options, default_val

# 2. Update Time Slider
@dash.callback(
    [Output('wave-time-slider', 'max'),
     Output('wave-time-slider', 'marks'),
     Output('wave-time-slider', 'value')],
    Input('wave-run-selector', 'value')
)
def update_time_slider(prefix):
    steps = get_available_timesteps(prefix)
    if not steps:
        return 0, {0: '0'}, 0
    marks = {s: str(s) for s in steps[::max(1, len(steps)//10)]}
    return max(steps), marks, min(steps)

# 3. Main Rendering Logic
@dash.callback(
    [Output('wave-tab-content', 'children'),
     Output('wave-z-slider', 'max'),
     Output('wave-z-slider', 'marks')],
    [Input('wave-tabs', 'active_tab'),
     Input('wave-run-selector', 'value'),
     Input('wave-time-slider', 'value'),
     Input('wave-var-selector', 'value'),
     Input('wave-z-slider', 'value')]
)
def render_wave_content(active_tab, prefix, t_step, var_name, z_layer):
    if not prefix:
        return html.Div("Please select a simulation project."), 0, {0: '0'}
    
    try:
        data = load_sim_snapshot(prefix, t_step)
        nz = data['grid']['nz']
        marks = {i: str(i) for i in range(nz)}
        
        # Color Scale: Diverging is best for waves
        color_scale = "RdBu" 

        if active_tab == "tab-2d":
            field_key = var_name
            if field_key not in data['fields']:
                field_key = list(data['fields'].keys())[0]
            
            field = data['fields'][field_key]
            slice_data = field[min(z_layer, nz-1), :, :]
            
            # Centralize around zero for diverging scale
            max_abs = float(max(abs(np.min(slice_data)), abs(np.max(slice_data))))
            if max_abs == 0: max_abs = 1.0

            fig = px.imshow(
                slice_data,
                color_continuous_scale=color_scale,
                zmin=-max_abs, zmax=max_abs,
                labels=dict(x="X", y="Y", color=var_name),
                title=f"{var_name} Field - Step {t_step} (Layer {z_layer})"
            )
            fig.update_layout(
                template="plotly_dark",
                paper_bgcolor="rgba(0,0,0,0)",
                plot_bgcolor="rgba(0,0,0,0)",
                margin=dict(l=20, r=20, t=40, b=20)
            )
            return dcc.Graph(figure=fig, style={"height": "600px"}), nz-1, marks
            
        elif active_tab == "tab-3d":
            # Add 3D Tab for consistency even if not in original layout, 
            # or skip if Wave layout has only 2 tabs. 
            # Wave layout currently has tab-2d and tab-trend.
            return html.Div("3D Volumetric Wavefield Analysis under development."), nz-1, marks
            
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
                template="plotly_dark",
                title=f"Evolution of {var_name} Field"
            )
            fig.add_vline(x=t_step, line_dash="dash", line_color="orange")
            fig.update_layout(
                paper_bgcolor="rgba(0,0,0,0)",
                plot_bgcolor="rgba(0,0,0,0)"
            )
            return dcc.Graph(figure=fig, style={"height": "600px"}), nz-1, marks
            
    except Exception as e:
        return html.Div(f"Extraction Error: {e}", className="text-danger p-4"), 0, {0: '0'}
