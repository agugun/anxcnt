from dash import dcc, html
import dash_bootstrap_components as dbc
from . import callbacks # Register callbacks

def layout():
    return dbc.Container([
        dbc.Row([
            dbc.Col([
                html.Div([
                    html.H2("Thermodynamics Analyzer", className="mb-1"),
                    html.P("Real-time monitoring of thermal diffusion fields.", className="text-muted"),
                ], className="mb-4")
            ], width=12)
        ]),
        
        dbc.Row([
            # Control Sidebar
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader([html.I(className="bi bi-sliders me-2"), "Simulation Controls"]),
                    dbc.CardBody([
                        html.Label("Target Project:", className="small text-muted mb-1"),
                        dcc.Dropdown(id='heat-run-selector', options=[], className="mb-4 custom-dropdown"),
                        
                        html.Label("Variable to Display:", className="small text-muted mb-1"),
                        dcc.Dropdown(
                            id='heat-var-selector',
                            options=['Temperature'],
                            value='Temperature',
                            className="mb-4 custom-dropdown"
                        ),
                        
                        html.Label("Temporal Step:", className="small text-muted mb-1"),
                        dcc.Slider(id='heat-time-slider', min=0, max=0, step=1, value=0, className="mb-4"),
                        
                        html.Label("Spatial Layer (Z):", className="small text-muted mb-1"),
                        dcc.Slider(id='heat-z-slider', min=0, max=0, step=1, value=0),
                        
                        dbc.Button(
                            [html.I(className="bi bi-arrow-clockwise me-2"), "Synchronize Data"], 
                            id='heat-btn-refresh', 
                            color="primary", 
                            className="w-100 mt-4"
                        )
                    ])
                ], className="h-100")
            ], width=3),
            
            # Module Main Render
            dbc.Col([
                dbc.Tabs([
                    dbc.Tab(label="Cross-Section", tab_id="tab-2d", label_class_name="fw-bold"),
                    dbc.Tab(label="Volumetric View", tab_id="tab-3d", label_class_name="fw-bold"),
                    dbc.Tab(label="Temporal Trends", tab_id="tab-trend", label_class_name="fw-bold")
                ], id="heat-tabs", active_tab="tab-2d", className="mb-0 custom-tabs"),
                
                html.Div(
                    id="heat-tab-content", 
                    className="p-4 bg-dark-soft border-glass-bottom", 
                    style={"height": "700px", "borderRadius": "0 0 16px 16px"}
                )
            ], width=9)
        ])
    ], fluid=True)
