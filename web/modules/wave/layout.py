from dash import dcc, html
import dash_bootstrap_components as dbc
from . import callbacks 

def layout():
    return dbc.Container([
        dbc.Row([
            dbc.Col([
                html.Div([
                    html.H2("Wave Propagation Module", className="mb-1"),
                    html.P("Dynamic visualization of wavefield evolution.", className="text-muted"),
                ], className="mb-4")
            ], width=12)
        ]),
        
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader([html.I(className="bi bi-broadcast me-2"), "Wavefield Controls"]),
                    dbc.CardBody([
                        html.Label("Project Selection:", className="small text-muted mb-1"),
                        dcc.Dropdown(id='wave-run-selector', options=[], className="mb-4 custom-dropdown"),
                        
                        html.Label("Field Variable:", className="small text-muted mb-1"),
                        dcc.Dropdown(
                            id='wave-var-selector',
                            options=['Displacement', 'Velocity'],
                            value='Displacement',
                            className="mb-4 custom-dropdown"
                        ),
                        
                        html.Label("Time Step:", className="small text-muted mb-1"),
                        dcc.Slider(id='wave-time-slider', min=0, max=0, step=1, value=0, className="mb-4"),
                        
                        html.Label("Layer (Z):", className="small text-muted mb-1"),
                        dcc.Slider(id='wave-z-slider', min=0, max=0, step=1, value=0),
                        
                        dbc.Button(
                            [html.I(className="bi bi-arrow-repeat me-2"), "Reload Data"], 
                            id='wave-btn-refresh', 
                            color="primary", 
                            className="w-100 mt-4"
                        )
                    ])
                ], className="h-100")
            ], width=3),
            
            dbc.Col([
                dbc.Tabs([
                    dbc.Tab(label="Snapshot", tab_id="tab-2d", label_class_name="fw-bold"),
                    dbc.Tab(label="Waveform History", tab_id="tab-trend", label_class_name="fw-bold")
                ], id="wave-tabs", active_tab="tab-2d", className="mb-0 custom-tabs"),
                
                html.Div(
                    id="wave-tab-content", 
                    className="p-4 bg-dark-soft border-glass-bottom", 
                    style={"height": "700px", "borderRadius": "0 0 16px 16px"}
                )
            ], width=9)
        ])
    ], fluid=True)
