import json
import numpy as np
import plotly.graph_objs as go
import plotly.offline as pyo

# ================================================================
# CONFIGURATION
# ================================================================
JSON_FILE = "output_results.json"     # your computed JSON file
OUTPUT_HTML = "coverage_3d_view.html" # generated 3D map


# ================================================================
# LOAD JSON
# ================================================================
print(f"\nLoading JSON: {JSON_FILE}")

with open(JSON_FILE, "r") as f:
    data = json.load(f)

if len(data) == 0:
    raise ValueError("JSON file is empty!")

print(f"Loaded {len(data)} samples.")


# ================================================================
# EXTRACT FIELDS
# ================================================================
lats = np.array([d["lat"] for d in data])
lons = np.array([d["lon"] for d in data])
elevs = np.array([d["terrain_elev_m"] for d in data])
clearance = np.array([d["clearance_m"] for d in data])
loss = np.array([d["knife_edge_loss_db"] for d in data])
dist = np.array([d["distance_m"] for d in data])


# ================================================================
# DETECT OBSTACLES
# ================================================================
obs_mask = (clearance < 0) | (loss > 0.5)

obs_lats = lats[obs_mask]
obs_lons = lons[obs_mask]
obs_elevs = elevs[obs_mask]

print(f"Obstacles detected: {len(obs_lats)}")


# ================================================================
# GUESS TOWER POSITION
# (first entry is always tower origin)
# ================================================================
tower_lat = lats[0]
tower_lon = lons[0]
tower_elev = elevs[0]

# ================================================================
# BUILD 3D OBJECTS
# ================================================================

# Terrain scatter
terrain_points = go.Scatter3d(
    x=lons,
    y=lats,
    z=elevs,
    mode='markers',
    marker=dict(
        size=2,
        color=elevs,
        colorscale='Viridis'
    ),
    name="Terrain"
)

# Obstacles
obstacle_points = go.Scatter3d(
    x=obs_lons,
    y=obs_lats,
    z=obs_elevs,
    mode='markers',
    marker=dict(
        size=4,
        color="red"
    ),
    name="Obstacles"
)

# Tower
tower_point = go.Scatter3d(
    x=[tower_lon],
    y=[tower_lat],
    z=[tower_elev + 10],  # slight elevation
    mode='markers',
    marker=dict(size=8, color="blue"),
    name="Tower"
)

# ================================================================
# LAYOUT
# ================================================================
layout = go.Layout(
    title="3D Map from JSON",
    scene=dict(
        xaxis_title="Longitude",
        yaxis_title="Latitude",
        zaxis_title="Elevation (m)",
        aspectmode='data'
    ),
    margin=dict(l=0, r=0, t=40, b=0)
)

# ================================================================
# FINAL FIGURE
# ================================================================
fig = go.Figure(data=[terrain_points, obstacle_points, tower_point], layout=layout)
pyo.plot(fig, filename=OUTPUT_HTML, auto_open=True)

print(f"\n✔ 3D Visualization saved as: {OUTPUT_HTML}\n")
