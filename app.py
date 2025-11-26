# app.py
import math
import json
import numpy as np
from pyproj import Geod
import srtm
from tqdm import tqdm
import plotly.graph_objs as go

import dash
from dash import dcc, html, Input, Output, State
import dash_bootstrap_components as dbc

# -----------------------------
# Constants for Earth curvature (effective radius)
# -----------------------------
R_EARTH = 6371000.0
K_FACTOR = 4.0 / 3.0
R_EFF = R_EARTH * K_FACTOR

# -----------------------------
# Knife-edge loss
# -----------------------------
def knife_edge_loss_db(v):
    if v <= -0.78:
        return 0.0
    return 6.9 + 20 * math.log10(math.sqrt((v - 0.1) ** 2 + 1) + v - 0.1)

def safe_elev(e):
    return float(e) if e is not None else 0.0

# -----------------------------
# Fresnel radius (using freq in Hz)
# -----------------------------
def fresnel_radius_freq(freq_hz, d1, d2):
    if d1 + d2 == 0:
        return 0.0
    lam = 299792458.0 / freq_hz
    return math.sqrt((lam * d1 * d2) / (d1 + d2))

# -----------------------------
# Earth bulge
# -----------------------------
def earth_bulge(d, D):
    return d * (D - d) / (2.0 * R_EFF)

# -----------------------------
# FSPL helper: compute max distance (km)
# -----------------------------
def compute_fspl_max_distance_km(freq_mhz, fspl_limit_db):
    return 10 ** ((fspl_limit_db - 20 * math.log10(freq_mhz) - 32.45) / 20)

# -----------------------------
# Core coverage computation (same logic you had)
# -----------------------------
def compute_coverage(lat, lon, height_m, freq_mhz,
                     az_step_deg=1, max_dist_km=10, sample_m=30,
                     rx_h=1.5, fresnel_frac=0.6, sensitivity_db=25):
    freq_hz = freq_mhz * 1e6
    geod = Geod(ellps="WGS84")

    max_dist = max_dist_km * 1000.0
    dists = np.arange(sample_m, max_dist + sample_m, sample_m)

    dem = srtm.get_data()
    tower_ground = safe_elev(dem.get_elevation(lat, lon))
    tower_abs = tower_ground + height_m

    results = []
    full_samples = []

    for az in tqdm(range(0, 360, int(az_step_deg)), desc="Azimuths"):
        max_usable_D = 0.0
        max_usable_mode = "LOS"
        max_usable_loss = 0.0
        max_usable_coord = (lat, lon)
        blocked = False

        for i, D in enumerate(dists):
            lonR, latR, _ = geod.fwd(lon, lat, az, D)
            elevR = safe_elev(dem.get_elevation(latR, lonR))
            rx_abs = elevR + rx_h

            D_usable = True
            D_max_loss = 0.0
            D_mode = "LOS"

            for j in range(i):
                d = dists[j]
                lonI, latI, _ = geod.fwd(lon, lat, az, d)
                elevI = safe_elev(dem.get_elevation(latI, lonI))

                line_h = tower_abs + (rx_abs - tower_abs) * (d / D)
                bulge = earth_bulge(d, D)

                d1 = d
                d2 = D - d
                if d2 <= 0 or (d1 + d2) == 0:
                    continue

                f1 = fresnel_radius_freq(freq_hz, d1, d2)
                clearance_required = f1 * fresnel_frac

                obst = elevI + bulge

                if obst >= (line_h - clearance_required):
                    h = obst - line_h
                    lam = 299792458.0 / freq_hz
                    denom = lam * d1 * d2
                    if denom == 0:
                        v = 0.0
                    else:
                        v = h * math.sqrt(2.0 * (d1 + d2) / denom)

                    loss = knife_edge_loss_db(v)
                    D_max_loss = max(D_max_loss, loss)

                    if loss > sensitivity_db:
                        D_usable = False
                        D_mode = "Blocked"
                        break
                    else:
                        D_mode = "Diffraction"

            full_samples.append({
                "azimuth": int(az),
                "distance_m": float(D),
                "lat": float(latR),
                "lon": float(lonR),
                "terrain_elev_m": float(elevR),
                "rx_abs_m": float(rx_abs),
                "usable": bool(D_usable),
                "mode": str(D_mode),
                "max_loss_db_on_path": float(D_max_loss)
            })

            if not D_usable:
                blocked = True
                break
            else:
                max_usable_D = D
                max_usable_mode = D_mode
                max_usable_loss = D_max_loss
                max_usable_coord = (latR, lonR)

        if not blocked and max_usable_D >= max_dist - 1e-6:
            max_usable_mode = "LOS"
            max_usable_loss = 0.0

        results.append({
            "azimuth": int(az),
            "max_distance_m": float(max_usable_D),
            "max_distance_km": float(max_usable_D / 1000.0),
            "coord_lat": float(max_usable_coord[0]),
            "coord_lon": float(max_usable_coord[1]),
            "mode": str(max_usable_mode),
            "loss_db": float(max_usable_loss)
        })

    return results, full_samples

# -----------------------------
# 3D Terrain builder
# -----------------------------
def build_3d_map(center_lat, center_lon, radius_km, elev, grid=200):
    GRID = grid
    lats = []
    lons = []
    heights = []

    deg_km = radius_km / 111.0

    lat_min = center_lat - deg_km
    lat_max = center_lat + deg_km
    lon_min = center_lon - deg_km
    lon_max = center_lon + deg_km

    for i in range(GRID):
        lat_row = []
        lon_row = []
        h_row = []
        for j in range(GRID):
            cur_lat = lat_min + (lat_max - lat_min) * (i / GRID)
            cur_lon = lon_min + (lon_max - lon_min) * (j / GRID)
            lat_row.append(cur_lat)
            lon_row.append(cur_lon)
            h_row.append(safe_elev(elev.get_elevation(cur_lat, cur_lon)))
        lats.append(lat_row)
        lons.append(lon_row)
        heights.append(h_row)

    return np.array(lats), np.array(lons), np.array(heights)

# -----------------------------
# Profiles (same defaults)
# -----------------------------
PROFILES = {
    "2G": {"freq_mhz": 900.0,  "fspl_limit_db": 120.0, "tower_height_m": 45.0, "rx_h": 1.5, "fresnel_frac": 0.6},
    "3G": {"freq_mhz": 2100.0, "fspl_limit_db": 115.0, "tower_height_m": 35.0, "rx_h": 1.5, "fresnel_frac": 0.6},
    "4G": {"freq_mhz": 1800.0, "fspl_limit_db": 125.0, "tower_height_m": 30.0, "rx_h": 1.5, "fresnel_frac": 0.6},
    "5G": {"freq_mhz": 3500.0, "fspl_limit_db": 110.0, "tower_height_m": 25.0, "rx_h": 1.5, "fresnel_frac": 0.5},
}

# -----------------------------
# Dash app layout
# -----------------------------
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server

app.layout = dbc.Container([
    dbc.Row([
        dbc.Col(html.H2("Tower Coverage — 3D Terrain Viewer"), width=8),
        dbc.Col(html.Div("Click a red point to see details →"), width=4, style={"textAlign":"right", "paddingTop":"10px"})
    ], align="center", className="my-2"),

    dbc.Row([
        dbc.Col([
            dbc.Card([
                dbc.CardBody([
                    dbc.Row([
                        dbc.Col([
                            dbc.Label("Latitude"),
                            dcc.Input(id="lat-input", type="number", value=18.5, step=0.000001, style={"width":"100%"}),
                        ], width=3),
                        dbc.Col([
                            dbc.Label("Longitude"),
                            dcc.Input(id="lon-input", type="number", value=80.9, step=0.000001, style={"width":"100%"}),
                        ], width=3),
                        dbc.Col([
                            dbc.Label("Network/Profile"),
                            dcc.Dropdown(id="profile-dropdown", options=[{"label":k,"value":k} for k in PROFILES.keys()], value="4G"),
                        ], width=2),
                        dbc.Col([
                            dbc.Label("Tower height (m)"),
                            dcc.Input(id="tower-height", type="number", value=30, step=1, style={"width":"100%"}),
                        ], width=2),
                        dbc.Col([
                            dbc.Label("Receiver height (m)"),
                            dcc.Input(id="rx-height", type="number", value=1.5, step=0.1, style={"width":"100%"}),
                        ], width=2),
                    ]),
                    html.Hr(),
                    dbc.Row([
                        dbc.Col([
                            dbc.Label("Azimuth step (deg)"),
                            dcc.Input(id="az-step", type="number", value=5, step=1, style={"width":"100%"}),
                        ], width=2),
                        dbc.Col([
                            dbc.Label("Sample (m)"),
                            dcc.Input(id="sample-m", type="number", value=30, step=1, style={"width":"100%"}),
                        ], width=2),
                        dbc.Col([
                            dbc.Label("Sensitivity (dB)"),
                            dcc.Input(id="sensitivity-db", type="number", value=25, step=1, style={"width":"100%"}),
                        ], width=2),
                        dbc.Col([
                            dbc.Label("Max distance (km) (auto if blank)"),
                            dcc.Input(id="max-dist-km", type="number", placeholder="auto", style={"width":"100%"}),
                        ], width=2),
                        dbc.Col([
                            dbc.Button("Run analysis", id="run-btn", color="primary", className="mt-2"),
                        ], width=2),
                    ]),
                    html.Div(id="run-status", className="mt-2")
                ])
            ])
        ], width=12),
    ], className="mb-3"),

    dbc.Row([
        dbc.Col([
            dcc.Loading(dcc.Graph(id="terrain-3d", figure={}, style={"height":"80vh"}), type="default"),
            dcc.Store(id="store-full-samples"),
            dcc.Store(id="store-results")
        ], width=9),

        dbc.Col([
            dbc.Card([
                dbc.CardBody([
                    html.H5("Point Info", className="card-title"),
                    html.Div(id="point-info", children=[
                        html.P("Click a blocking point on the 3D map to see details here.")
                    ])
                ])
            ])
        ], width=3)
    ], className="mb-3"),
], fluid=True)

# -----------------------------
# Callbacks
# -----------------------------
@app.callback(
    Output("run-status", "children"),
    Output("terrain-3d", "figure"),
    Output("store-full-samples", "data"),
    Output("store-results", "data"),
    Input("run-btn", "n_clicks"),
    State("lat-input", "value"),
    State("lon-input", "value"),
    State("profile-dropdown", "value"),
    State("tower-height", "value"),
    State("rx-height", "value"),
    State("az-step", "value"),
    State("sample-m", "value"),
    State("sensitivity-db", "value"),
    State("max-dist-km", "value"),
    prevent_initial_call=True
)
def run_analysis(n_clicks, lat, lon, profile_key, tower_height, rx_height, az_step, sample_m, sensitivity_db, max_dist_km):
    try:
        status = "Loading SRTM data..."
        # prepare profile defaults
        prof = PROFILES.get(profile_key, PROFILES["4G"])
        freq_mhz = prof["freq_mhz"]
        fspl_default = prof["fspl_limit_db"]

        # compute auto max distance if not set
        if max_dist_km is None or max_dist_km == "" or float(max_dist_km) <= 0:
            max_dist_km_use = compute_fspl_max_distance_km(freq_mhz, fspl_default)
            if max_dist_km_use < 0.1:
                max_dist_km_use = 0.1
        else:
            max_dist_km_use = float(max_dist_km)

        status = f"Computing coverage (max_dist={max_dist_km_use:.3f} km, freq={freq_mhz} MHz)..."

        elev = srtm.get_data()

        # run heavy compute
                # sanitize inputs (avoid None)
        lat = float(lat or 0.0)
        lon = float(lon or 0.0)
        tower_height = float(tower_height or prof["tower_height_m"])
        rx_height = float(rx_height or prof["rx_h"])
        az_step = int(az_step or 1)
        sample_m = float(sample_m or 30.0)
        sensitivity_db = float(sensitivity_db or 25.0)

        results, full_samples = compute_coverage(
            lat, lon, tower_height, float(freq_mhz),
            az_step_deg=az_step, max_dist_km=float(max_dist_km_use),
            sample_m=sample_m, rx_h=rx_height,
            fresnel_frac=float(prof["fresnel_frac"]),
            sensitivity_db=sensitivity_db
        )


        # Save outputs
        with open("coverage_direction_report.json", "w") as f:
            json.dump(results, f, indent=2)
        with open("output_results.json", "w") as f:
            json.dump(full_samples, f, indent=2)

        # Build 3D surface
        lats, lons, heights = build_3d_map(float(lat), float(lon), float(max_dist_km_use), elev, grid=200)

        surface = go.Surface(
            z=heights,
            x=lons,
            y=lats,
            showscale=False,
            opacity=0.9,
            name="Terrain"
        )

        # gather blocked points
        blocked_points = [s for s in full_samples if s.get("usable") is False]
        obs_lat = [float(s["lat"]) for s in blocked_points]
        obs_lon = [float(s["lon"]) for s in blocked_points]
        obs_hgt = [float(s["terrain_elev_m"]) for s in blocked_points]
        obs_dist_km = [float(s["distance_m"]) / 1000.0 for s in blocked_points]
        obs_az = [int(s["azimuth"]) for s in blocked_points]
        obs_mode = [s.get("mode", "") for s in blocked_points]
        obs_loss = [float(s.get("max_loss_db_on_path", 0.0)) for s in blocked_points]

        obstacle_points = go.Scatter3d(
            x=obs_lon,
            y=obs_lat,
            z=obs_hgt,
            mode='markers',
            marker=dict(size=4, color="red"),
            name="Blocking Points",
            customdata=np.stack([obs_lat, obs_lon, obs_dist_km, obs_az, obs_mode, obs_loss], axis=1).tolist(),
            hovertemplate="<b>Blocked</b><br>Lat: %{customdata[0]:.6f}<br>Lon: %{customdata[1]:.6f}<br>Dist: %{customdata[2]:.3f} km<br>Az: %{customdata[3]}°<extra></extra>"
        )

        tower_ground = safe_elev(elev.get_elevation(float(lat), float(lon)))
        tower_marker = go.Scatter3d(
            x=[float(lon)],
            y=[float(lat)],
            z=[tower_ground + float(tower_height)],
            mode='markers',
            marker=dict(size=8, color="blue"),
            name="Tower",
            customdata=[[float(lat), float(lon), 0.0, -1, "Tower", 0.0]],
            hovertemplate="<b>Tower</b><br>Lat: %{customdata[0]:.6f}<br>Lon: %{customdata[1]:.6f}<extra></extra>"
        )

        fig = go.Figure(data=[surface, obstacle_points, tower_marker])
        fig.update_layout(
            title=f"3D Terrain Map with Blocking Points (center: {lat:.6f},{lon:.6f})",
            scene=dict(
                xaxis_title="Longitude",
                yaxis_title="Latitude",
                zaxis_title="Elevation (m)"
            ),
            margin=dict(l=0, r=0, b=0, t=40),
            annotations=[dict(x=0.5, y=1.13, xref='paper', yref='paper', showarrow=False,
                              text="N ↑ | NE ↗ | E → | SE ↘ | S ↓ | SW ↙ | W ← | NW ↖",
                              font=dict(size=13))]
        )

        status = f"Done — {len(blocked_points)} blocking points found. Click a red point to see details."

        # store full_samples and results as JSON for client-side callbacks
        return status, fig, full_samples, results

    except Exception as e:
        return f"Error: {str(e)}", dash.no_update, dash.no_update, dash.no_update

# Update point info when user clicks a point in the graph
@app.callback(
    Output("point-info", "children"),
    Input("terrain-3d", "clickData"),
    State("store-full-samples", "data"),
    prevent_initial_call=False
)
def display_point_info(clickData, full_samples):
    # if no clickData, show default message
    if not clickData:
        return html.P("Click a blocking point on the 3D map to see details here.")

    try:
        pt = clickData["points"][0]
        custom = pt.get("customdata", None)
        if not custom:
            # Possibly clicked on surface; show coordinates
            x = pt.get("x", None)
            y = pt.get("y", None)
            z = pt.get("z", None)
            return html.Div([
                html.P(f"Clicked on surface / non-blocking point"),
                html.P(f"Longitude: {x}"),
                html.P(f"Latitude: {y}"),
                html.P(f"Elevation: {z}")
            ])

        lat = float(custom[0])
        lon = float(custom[1])
        dist_km = float(custom[2])
        az = int(custom[3])
        mode = str(custom[4]) if len(custom) >= 5 else ""
        loss = float(custom[5]) if len(custom) >= 6 else 0.0

        return html.Div([
            html.P([html.B("Latitude: "), f"{lat:.6f}"]),
            html.P([html.B("Longitude: "), f"{lon:.6f}"]),
            html.P([html.B("Distance (km): "), f"{dist_km:.3f}"]),
            html.P([html.B("Azimuth (°): "), f"{az}"]),
            html.P([html.B("Mode: "), f"{mode}"]),
            html.P([html.B("Max loss (dB): "), f"{loss:.2f}"]),
        ])
    except Exception as e:
        return html.P(f"Error reading clicked point: {str(e)}")


if __name__ == "__main__":
    app.run(debug=True)
