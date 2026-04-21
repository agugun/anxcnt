from dash import dcc, html
import dash_bootstrap_components as dbc
from . import callbacks 

def layout():
    return dbc.Container([
        dbc.Row([
            dbc.Col([
                html.Div([
                    html.H2("Pressure Analysis Module", className="mb-1"),
                    html.P("Geomechanical and fluid pressure equilibrium monitoring.", className="text-muted"),
                ], className="mb-4")
            ], width=12)
        ]),
        
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader([html.I(className="bi bi-speedometer2 me-2"), "Pressure Controls"]),
                    dbc.CardBody([
                        html.Label("Project Selection:", className="small text-muted mb-1"),
                        dcc.Dropdown(id='pressure-run-selector', options=[], className="mb-4 custom-dropdown"),
                        
                        html.Label("Variable:", className="small text-muted mb-1"),
                        dcc.Dropdown(
                            id='pressure-var-selector',
                            options=['Pressure'],
                            value='Pressure',
                            className="mb-4 custom-dropdown"
                        ),
                        
                        html.Label("Time Step:", className="small text-muted mb-1"),
                        dcc.Slider(id='pressure-time-slider', min=0, max=0, step=1, value=0, className="mb-4"),
                        
                        html.Label("Layer (Z):", className="small text-muted mb-1"),
                        dcc.Slider(id='pressure-z-slider', min=0, max=0, step=1, value=0),
                        
                        dbc.Button(
                            [html.I(className="bi bi-arrow-repeat me-2"), "Refresh Data"], 
                            id='pressure-btn-refresh', 
                            color="primary", 
                            className="w-100 mt-4"
                        )
                    ])
                ], className="h-100")
            ], width=3),
            
            dbc.Col([
                dbc.Tabs([
                    dbc.Tab(label="2D Gradient", tab_id="tab-2d", label_class_name="fw-bold"),
                    dbc.Tab(label="3D Matrix", tab_id="tab-3d", label_class_name="fw-bold"),
                    dbc.Tab(label="History Plot", tab_id="tab-trend", label_class_name="fw-bold")
                ], id="pressure-tabs", active_tab="tab-2d", className="mb-0 custom-tabs"),
                
                html.Div(
                    id="pressure-tab-content", 
                    className="p-4 bg-dark-soft border-glass-bottom", 
                    style={"height": "700px", "borderRadius": "0 0 16px 16px"}
                )
            ], width=9)
        ])
    ], fluid=True)
