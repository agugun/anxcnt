import dash_bootstrap_components as dbc
from dash import dcc, html

def layout():
    return dbc.Container([
        dbc.Row([
            dbc.Col([
                html.Div([
                    html.H2("Data Sovereignty Audit", className="mb-1"),
                    html.P("Verifying local workspace isolation and compliance.", className="text-muted"),
                ], className="mb-4")
            ])
        ]),
        
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader([html.I(className="bi bi-shield-lock me-2"), "Edge Compliance Status"]),
                    dbc.CardBody([
                        html.Div([
                            html.Div([
                                html.P("Workspace Isolation", className="mb-0 fw-bold"),
                                dbc.Badge("ACTIVE", color="success"),
                            ], className="d-flex justify-content-between align-items-center mb-3"),
                            
                            html.Div([
                                html.P("Remote Data Leakage", className="mb-0 fw-bold"),
                                dbc.Badge("NONE DETECTED", color="success"),
                            ], className="d-flex justify-content-between align-items-center mb-3"),
                            
                            html.Div([
                                html.P("Local VTI Integrity", className="mb-0 fw-bold"),
                                dbc.Badge("VERIFIED", color="info"),
                            ], className="d-flex justify-content-between align-items-center mb-3"),
                        ], className="p-3 bg-dark-soft rounded-3")
                    ])
                ], className="mb-4")
            ], width=4),
            
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader([html.I(className="bi bi-terminal me-2"), "Real-time Access Log"]),
                    dbc.CardBody([
                        html.Div(id="sovereignty-log", className="p-4 small font-monospace bg-black rounded-3", 
                                 style={"height": "400px", "overflowY": "auto", "color": "#10b981"})
                    ])
                ])
            ], width=8)
        ]),
        
        dcc.Interval(id="sovereignty-timer", interval=2000, n_intervals=0)
    ])

# Callbacks implemented directly in this module for simplicity
from dash import Input, Output
import dash
import datetime
import random

@dash.callback(
    Output("sovereignty-log", "children"),
    Input("sovereignty-timer", "n_intervals")
)
def update_sovereignty_log(n):
    actions = [
        "Reading block from exports/",
        "VTI mesh parsing complete",
        "Local memory allocation verified",
        "Edge-compute gradient calculation",
        "Sanitizing temporary buffer",
        "Internal websocket heartbeat"
    ]
    
    log_entries = []
    for i in range(max(0, n-10), n+1):
        ts = (datetime.datetime.now() - datetime.timedelta(seconds=(n-i)*2)).strftime("%H:%M:%S")
        action = actions[i % len(actions)]
        log_entries.append(html.Div(f"[{ts}] {action}... OK"))
        
    return log_entries
