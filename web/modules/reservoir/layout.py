from dash import dcc, html
import dash_bootstrap_components as dbc

# Import callbacks to register them
from . import callbacks

def layout():
    return dbc.Container([
        dbc.Row([
            dbc.Col([
                html.Div([
                    html.H2("Reservoir Monitoring Module", className="mb-1"),
                    html.P("Real-time surveillance of subsurface dynamics and volumetric trends.", className="text-muted"),
                ], className="mb-4")
            ], width=12)
        ]),
        
        dbc.Row([
            # Control Sidebar
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader([html.I(className="bi bi-gear-wide-connected me-2"), "Reservoir Controls"]),
                    dbc.CardBody([
                        html.Label("Simulation Case:", className="small text-muted mb-1"),
                        dcc.Dropdown(id='run-selector', options=[], className="mb-4 custom-dropdown"),
                        
                        html.Label("Variable:", className="small text-muted mb-1"),
                        dcc.Dropdown(
                            id='var-selector',
                            options=['Pressure', 'Sw', 'Sg'],
                            value='Pressure',
                            className="mb-4 custom-dropdown"
                        ),
                        
                        html.Label("Time Step:", className="small text-muted mb-1"),
                        dcc.Slider(id='time-slider', min=0, max=0, step=1, value=0, className="mb-4"),
                        
                        html.Label("Layer (Z):", className="small text-muted mb-1"),
                        dcc.Slider(id='z-slider', min=0, max=0, step=1, value=0),
                        
                        dbc.Button(
                            [html.I(className="bi bi-arrow-clockwise me-2"), "Refresh Runs"], 
                            id='btn-refresh', 
                            color="primary", 
                            className="w-100 mt-4"
                        )
                    ])
                ], className="h-100")
            ], width=3),
            
            # Module Main View
            dbc.Col([
                dbc.Tabs([
                    dbc.Tab(label="2D Slice", tab_id="tab-2d", label_class_name="fw-bold"),
                    dbc.Tab(label="3D Reservoir", tab_id="tab-3d", label_class_name="fw-bold"),
                    dbc.Tab(label="Volumetric Trend", tab_id="tab-trend", label_class_name="fw-bold")
                ], id="tabs", active_tab="tab-2d", className="mb-0 custom-tabs"),
                
                html.Div(
                    id="tab-content", 
                    className="p-4 bg-dark-soft border-glass-bottom", 
                    style={"height": "700px", "borderRadius": "0 0 16px 16px"}
                )
            ], width=9)
        ])
    ], fluid=True)
