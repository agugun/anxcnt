import dash
from dash import dcc, html, Input, Output, State
import dash_bootstrap_components as dbc
import os
import sys

# Add project root to path for dynamic imports
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from registry import ModuleRegistry

# Initialize App
app = dash.Dash(__name__, 
                external_stylesheets=[dbc.themes.CYBORG, dbc.icons.BOOTSTRAP],
                suppress_callback_exceptions=True)
server = app.server

# Registry Setup
registry = ModuleRegistry()
registry.discover_modules()

# App Layout
app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    
    # Global Header
    dbc.Navbar(
        dbc.Container([
            html.A(
                dbc.Row([
                    dbc.Col(html.I(className="bi bi-intersect", style={"fontSize": "1.5rem", "color": "#6366f1"})),
                    dbc.Col(dbc.NavbarBrand("AXSCNT EDGE", className="ms-2 fw-bold", style={"letterSpacing": "0.1em"})),
                ], align="center", className="g-0"),
                href="/",
                style={"textDecoration": "none"},
            ),
            dbc.Button(
                [html.I(className="bi bi-folder2-open me-2"), "Local Workspace"],
                id="btn-connect-local",
                color="primary",
                size="sm",
                className="ms-auto"
            )
        ], fluid=True),
        color="rgba(10, 10, 12, 0.8)",
        dark=True,
        sticky="top",
        style={"backdropFilter": "blur(10px)", "borderBottom": "1px solid rgba(255,255,255,0.1)"}
    ),
    
    html.Div([
        # Dynamic Sidebar
        html.Div([
            html.Div([
                html.Small("ANALYTICS ENGINE", className="text-muted fw-bold px-4"),
            ], className="mt-4 mb-2"),
            dbc.Nav(
                [
                    dbc.NavLink(
                        [html.I(className="bi bi-house-door me-3"), "Dashboard"], 
                        href="/", active="exact"
                    )
                ] + [
                    dbc.NavLink(
                        [html.I(className="bi bi-activity me-3"), mod["label"]], 
                        href=f"/module/{mod['value']}", 
                        active="partial"
                    )
                    for mod in registry.get_module_list()
                ],
                vertical=True,
                pills=True,
                className="px-2"
            ),
            
            html.Div(style={"flexGrow": 1}),
            
            # Bottom section of sidebar
            html.Div([
                html.Hr(style={"backgroundColor": "rgba(255,255,255,0.1)"}),
                html.Div([
                    html.Small("v2.0.0-standardized", className="text-muted"),
                ], className="px-4 pb-4")
            ])
        ], className="sidebar-glass", style={
            "width": "260px", 
            "position": "fixed", 
            "top": "56px", 
            "bottom": 0, 
            "display": "flex",
            "flexDirection": "column"
        }),
        
        # Content Area
        html.Div(id="page-content", style={
            "marginLeft": "260px", 
            "padding": "2.5rem",
            "minHeight": "calc(100vh - 56px)"
        })
    ])
], style={"overflowX": "hidden"})

# Routing Callback
@app.callback(
    Output('page-content', 'children'),
    Input('url', 'pathname')
)
def display_page(pathname):
    if pathname == "/" or pathname is None:
        module_list = registry.get_module_list()
        return html.Div([
            html.Div([
                html.H1("Axscnt Edge Analytics", className="display-4 mb-2"),
                html.P("Empowering localized high-performance numerical intelligence.", className="lead text-muted"),
            ], className="mb-5"),
            
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardHeader([html.I(className="bi bi-cpu me-2"), "System Status"]),
                        dbc.CardBody([
                            html.H5(f"{len(registry.modules)} Modules Active", className="text-info"),
                            html.P("Local Workspace: /home/gugun/repo/axscnt", className="small text-muted"),
                            dbc.Badge("Engine Standardized", color="success", className="me-2"),
                            dbc.Badge("VTI Loader Active", color="primary"),
                        ])
                    ], className="h-100")
                ], width=4),
                
                dbc.Col([
                    dbc.Card([
                        dbc.CardHeader([html.I(className="bi bi-shield-check me-2"), "Data Sovereignty"]),
                        dbc.CardBody([
                            html.P("Computation and visualization are performing strictly on the edge. No simulation data leaves this workspace.", className="card-text"),
                            html.Div([
                                html.I(className="bi bi-arrow-down-short me-1"), "Latency: 0ms (Direct Disk Access)"
                            ], className="small text-success fw-bold")
                        ])
                    ], className="h-100")
                ], width=8),
            ], className="mb-5 g-4"),

            html.H4("Module Quick-Access", className="mb-4 mt-2 px-1 text-muted small fw-bold uppercase"),
            dbc.Row([
                dbc.Col([
                    html.A([
                        dbc.Card([
                            dbc.CardBody([
                                html.I(className=f"bi {mod.get('icon', 'bi-box')} mb-3 display-6", style={"color": "#6366f1"}),
                                html.H5(mod["label"], className="card-title"),
                                html.P(mod.get("description", "Launch domain analyzer"), className="card-text small text-muted"),
                            ], className="text-center p-4")
                        ], className="module-card-hover")
                    ], href=f"/module/{mod['value']}", style={"textDecoration": "none"})
                ], width=3) for mod in registry.get_module_list()
            ], className="g-4")
        ])
    
    if pathname.startswith("/module/"):
        module_id = pathname.split("/")[-1]
        mod_info = registry.modules.get(module_id)
        if mod_info and mod_info["layout"]:
            # If layout is a function, call it
            layout = mod_info["layout"]
            if callable(layout):
                return layout()
            return layout
            
    return html.Div([
        html.H1("404 - Not Found", className="text-danger"),
        html.P(f"The module at {pathname} was not found or is misconfigured.")
    ], className="text-center mt-5")

if __name__ == "__main__":
    app.run(debug=True, port=8050)
