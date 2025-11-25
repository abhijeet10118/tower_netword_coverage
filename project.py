# # """
# # tower_coverage.py
# # Takes user input for tower coordinates + height + frequency
# # Computes obstacle clearance using SRTM elevation data
# # """

# # import math
# # import numpy as np
# # from pyproj import Geod
# # import srtm
# # from tqdm import tqdm

# # # Constants
# # R_EARTH = 6371000.0
# # K_FACTOR = 4/3
# # R_EFF = R_EARTH * K_FACTOR


# # # Knife Edge Loss
# # def knife_edge_loss_db(v):
# #     if v <= -0.78:
# #         return 0.0
# #     return 6.9 + 20 * math.log10(math.sqrt((v - 0.1)**2 + 1) + v - 0.1)


# # def fresnel_radius(freq_hz, d1, d2):
# #     lam = 299792458.0 / freq_hz
# #     if d1 + d2 == 0:
# #         return 0
# #     return math.sqrt(lam * d1 * d2 / (d1 + d2))


# # def earth_bulge(d, D):
# #     return d * (D - d) / (2 * R_EFF)


# # # Load SRTM DEM
# # dem = srtm.get_data()


# # def get_elev(lat, lon):
# #     e = dem.get_elevation(lat, lon)
# #     return e if e is not None else 0.0


# # # MAIN COMPUTATION
# # def compute_coverage(lat, lon, height_m, freq_mhz,
# #                      az_step_deg=10, max_dist_km=10, sample_m=50):

# #     freq_hz = freq_mhz * 1e6
# #     geod = Geod(ellps="WGS84")

# #     max_dist = max_dist_km * 1000
# #     dists = np.arange(sample_m, max_dist, sample_m)

# #     tower_ground = get_elev(lat, lon)
# #     tower_abs = tower_ground + height_m

# #     results = []

# #     print("\nComputing coverage...")

# #     for az in tqdm(range(0, 360, az_step_deg)):
# #         blocked = None
# #         block_point = None
# #         block_loss = 0

# #         for i, D in enumerate(dists):
# #             lonR, latR, _ = geod.fwd(lon, lat, az, D)
# #             elevR = get_elev(latR, lonR)
# #             rx_abs = elevR + 1.5

# #             # Check all intermediate points
# #             for j in range(i):
# #                 d = dists[j]
# #                 lonI, latI, _ = geod.fwd(lon, lat, az, d)
# #                 elevI = get_elev(latI, lonI)

# #                 line_h = tower_abs + (rx_abs - tower_abs) * (d / D)
# #                 bulge = earth_bulge(d, D)

# #                 d1, d2 = d, D - d
# #                 if d2 <= 0:
# #                     continue

# #                 f1 = fresnel_radius(freq_hz, d1, d2)
# #                 clearance = f1 * 0.6
# #                 obst = elevI + bulge

# #                 if obst >= (line_h - clearance):
# #                     h = obst - line_h
# #                     lam = 299792458 / freq_hz
# #                     v = h * math.sqrt(2*(d1+d2)/(lam*d1*d2))
# #                     loss = knife_edge_loss_db(v)

# #                     mode = "Diffraction" if loss <= 20 else "Blocked"
# #                     blocked = (D, mode)
# #                     block_point = (latI, lonI)
# #                     block_loss = loss
# #                     break

# #             if blocked:
# #                 break

# #         if blocked:
# #             dist, mode = blocked
# #         else:
# #             dist, mode = max_dist, "LOS"

# #         results.append({
# #             "azimuth": az,
# #             "distance_m": dist,
# #             "mode": mode,
# #             "block_point": block_point,
# #             "loss_db": block_loss
# #         })

# #     return results


# # # ---------------- USER INPUT --------------------
# # if __name__ == "__main__":
# #     print("=== Cell Tower Obstacle & Coverage Tool ===")

# #     tower_lat = float(input("Enter tower latitude: "))
# #     tower_lon = float(input("Enter tower longitude: "))
# #     tower_height = float(input("Enter tower height above ground (m): "))
# #     freq_mhz = float(input("Enter frequency (MHz): "))

# #     output = compute_coverage(tower_lat, tower_lon, tower_height, freq_mhz)

# #     print("\n=== RESULTS ===")
# #     for r in output:
# #         km = r["distance_m"] / 1000
# #         print(f'Az {r["azimuth"]}° → {km:.2f} km → {r["mode"]}')
# #         if r["block_point"]:
# #             print(f'   Block at {r["block_point"]} | Loss: {r["loss_db"]:.1f} dB')
# import math
# import json
# import csv
# import numpy as np
# from pyproj import Geod
# import srtm
# from tqdm import tqdm
# import plotly.graph_objs as go
# import plotly.offline as pyo

# # -----------------------------
# # Knife-edge loss calculation
# # -----------------------------
# def knife_edge_loss_db(v):
#     if v <= -0.78:
#         return 0.0
#     return 6.9 + 20 * math.log10(math.sqrt((v - 0.1)**2 + 1) + v - 0.1)

# def safe_elev(e):
#     return e if e is not None else 0.0



# # -----------------------------
# # Fresnel Zone Radius
# # -----------------------------
# def fresnel_radius(d1, d2, freq_mhz):
#     f_hz = freq_mhz * 1e6
#     return math.sqrt((3e8 * d1 * d2) / (f_hz * (d1 + d2)))

# # -----------------------------
# # LOS / obstacle checking for each azimuth
# # -----------------------------
# def analyze_azimuth(lat, lon, tower_h, freq_mhz, azimuth_deg,
#                     max_dist_km, sample_m, rx_h, fresnel_frac, geod, elev):
    
#     az_results = []
#     dist_m = 0

#     while dist_m <= max_dist_km * 1000:

#         lon2, lat2, _ = geod.fwd(lon, lat, azimuth_deg, dist_m)

#         tx_ground = safe_elev(elev.get_elevation(lat, lon))
#         rx_ground = safe_elev(elev.get_elevation(lat2, lon2))

#         tx_height = tower_h + tx_ground
#         rx_height = rx_h + rx_ground

#         d1 = dist_m
#         d2 = 1
#         fresnel_r = fresnel_radius(d1, d2, freq_mhz) * fresnel_frac

#         los_line = tx_height - (tx_height - rx_height) * (dist_m / (max_dist_km * 1000))
#         terrain_h = rx_ground
#         clearance = terrain_h - los_line

#         v = (math.sqrt(2) * clearance) / fresnel_r if fresnel_r != 0 else 0
#         loss = knife_edge_loss_db(v)

#         az_results.append({
#             "lat": lat2,
#             "lon": lon2,
#             "distance_m": dist_m,
#             "terrain_elev_m": terrain_h,
#             "clearance_m": clearance,
#             "knife_edge_loss_db": loss
#         })

#         dist_m += sample_m

#     return az_results


# # -----------------------------
# # Build 3D elevation model using SRTM
# # -----------------------------


# def build_3d_map(center_lat, center_lon, radius_km, elev):
#     e = elev.get_elevation(lat, lon)
#     height = safe_elev(e)

#     GRID = 200  # resolution
#     lats = []
#     lons = []
#     heights = []

#     lat_min = center_lat - (radius_km / 111)
#     lat_max = center_lat + (radius_km / 111)
#     lon_min = center_lon - (radius_km / 111)
#     lon_max = center_lon + (radius_km / 111)

#     for i in range(GRID):
#         row_lat = []
#         row_lon = []
#         row_hgt = []
#         for j in range(GRID):
#             lat = lat_min + (lat_max - lat_min) * (i / GRID)
#             lon = lon_min + (lon_max - lon_min) * (j / GRID)

#             row_lat.append(lat)
#             row_lon.append(lon)
#             row_hgt.append(safe_elev(elev.get_elevation(lat, lon)))


#         lats.append(row_lat)
#         lons.append(row_lon)
#         heights.append(row_hgt)

#     return np.array(lats), np.array(lons), np.array(heights)


# # =================================================================
# # MAIN PROGRAM (TAKES NORMAL INPUT)
# # =================================================================

# print("\n=== Tower Obstacle Detection + 3D Map ===\n")

# lat = float(input("Enter tower latitude  : "))
# lon = float(input("Enter tower longitude : "))
# tower_height = float(input("Enter tower height (m) : "))
# freq_mhz = float(input("Enter frequency (MHz)  : "))
# az_step = float(input("Azimuth step (default 5°): ") or 5)
# max_dist_km = float(input("Max distance km (default 10 km): ") or 10)
# sample_m = float(input("Sampling distance meters (default 30m): ") or 30)
# rx_height = float(input("Receiver height (default 1.5m): ") or 1.5)
# fresnel_fraction = float(input("Fresnel fraction (default 0.6): ") or 0.6)

# print("\nLoading SRTM elevation data...")
# elev = srtm.get_data()

# geod = Geod(ellps="WGS84")
# all_samples = []
# obstacles = []

# print("\nProcessing azimuths...\n")

# for az in tqdm(range(0, 360, int(az_step))):
#     samples = analyze_azimuth(
#         lat, lon,
#         tower_height,
#         freq_mhz,
#         az,
#         max_dist_km,
#         sample_m,
#         rx_height,
#         fresnel_fraction,
#         geod,
#         elev
#     )
#     all_samples.extend(samples)

#     # detect obstacles
#     for s in samples:
#         if s["clearance_m"] < 0 or s["knife_edge_loss_db"] > 0:
#             obstacles.append(s)

# # Save JSON
# with open("output_results.json", "w") as f:
#     json.dump(all_samples, f, indent=2)

# print("\n✔ Saved calculations as output_results.json")

# # -----------------------------
# # 3D MAP GENERATION
# # -----------------------------
# print("\nGenerating 3D terrain model (this can take a minute)...")

# lats, lons, heights = build_3d_map(lat, lon, max_dist_km, elev)

# # 3D Surface
# surface = go.Surface(
#     z=heights,
#     x=lons,
#     y=lats,
#     showscale=False,
#     opacity=0.9
# )

# # Mark obstacles
# obs_lat = [o["lat"] for o in obstacles]
# obs_lon = [o["lon"] for o in obstacles]
# obs_hgt = [o["terrain_elev_m"] for o in obstacles]

# obstacle_points = go.Scatter3d(
#     x=obs_lon,
#     y=obs_lat,
#     z=obs_hgt,
#     mode='markers',
#     marker=dict(size=4, color="red"),
#     name="Obstacles"
# )

# # Tower marker
# # Tower marker (SAFE elevation)
# tower_ground = safe_elev(elev.get_elevation(lat, lon))

# tower_marker = go.Scatter3d(
#     x=[lon],
#     y=[lat],
#     z=[tower_ground + tower_height],    # FIXED: no None crash
#     mode='markers',
#     marker=dict(size=8, color="blue"),
#     name="Tower"
# )


# layout = go.Layout(
#     title="3D Terrain Map with Obstacles",
#     scene=dict(
#         xaxis_title="Longitude",
#         yaxis_title="Latitude",
#         zaxis_title="Elevation (m)"
#     ),
#     margin=dict(l=0, r=0, b=0, t=30)
# )

# fig = go.Figure(data=[surface, obstacle_points, tower_marker], layout=layout)

# output_html = "3d_map.html"
# pyo.plot(fig, filename=output_html, auto_open=False)

# print(f"\n✔ 3D interactive map saved as: {output_html}")
# print("Open it in any browser (Chrome/Edge) to rotate and view obstacles.\n")

# print("\n=== DONE ===\n")
import math
import json
import numpy as np
from pyproj import Geod
import srtm
from tqdm import tqdm
import plotly.graph_objs as go
import plotly.offline as pyo


# -----------------------------
# Knife-edge loss
# -----------------------------
def knife_edge_loss_db(v):
    if v <= -0.78:
        return 0.0
    return 6.9 + 20 * math.log10(math.sqrt((v - 0.1)**2 + 1) + v - 0.1)


def safe_elev(e):
    return e if e is not None else 0.0


# -----------------------------
# Fresnel Radius
# -----------------------------
def fresnel_radius(d1, d2, freq_mhz):
    f_hz = freq_mhz * 1e6
    return math.sqrt((3e8 * d1 * d2) / (f_hz * (d1 + d2)))


# -----------------------------
# LOS / obstacle analysis per azimuth
# -----------------------------
def analyze_azimuth(lat, lon, tower_h, freq_mhz, azimuth_deg,
                    max_dist_km, sample_m, rx_h, fresnel_frac, geod, elev):

    az_results = []
    dist_m = 0

    while dist_m <= max_dist_km * 1000:

        lon2, lat2, _ = geod.fwd(lon, lat, azimuth_deg, dist_m)

        tx_ground = safe_elev(elev.get_elevation(lat, lon))
        rx_ground = safe_elev(elev.get_elevation(lat2, lon2))

        tx_height = tower_h + tx_ground
        rx_height = rx_h + rx_ground

        d1 = dist_m
        d2 = 1
        fresnel_r = fresnel_radius(d1, d2, freq_mhz) * fresnel_frac

        if max_dist_km == 0:
            los_line = tx_height
        else:
            los_line = tx_height - (tx_height - rx_height) * (dist_m / (max_dist_km * 1000))

        terrain_h = rx_ground
        clearance = terrain_h - los_line

        v = (math.sqrt(2) * clearance) / fresnel_r if fresnel_r != 0 else 0
        loss = knife_edge_loss_db(v)

        az_results.append({
            "lat": lat2,
            "lon": lon2,
            "distance_m": dist_m,
            "terrain_elev_m": terrain_h,
            "clearance_m": clearance,
            "knife_edge_loss_db": loss
        })

        dist_m += sample_m

    return az_results


# -----------------------------
# FIXED 3D Terrain builder
# -----------------------------
def build_3d_map(center_lat, center_lon, radius_km, elev):

    GRID = 200  # high resolution
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


# =================================================================
# MAIN PROGRAM
# =================================================================
print("\n=== Tower Obstacle Detection + 3D Map ===\n")

lat = float(input("Enter tower latitude  : "))
lon = float(input("Enter tower longitude : "))
tower_height = float(input("Enter tower height (m) : "))
freq_mhz = float(input("Enter frequency (MHz)  : "))
az_step = float(input("Azimuth step (default 5°): ") or 5)
max_dist_km = float(input("Max distance km (default 10 km): ") or 10)
sample_m = float(input("Sampling distance meters (default 30m): ") or 30)
rx_height = float(input("Receiver height (default 1.5m): ") or 1.5)
fresnel_fraction = float(input("Fresnel fraction (default 0.6): ") or 0.6)

print("\nLoading SRTM elevation data...")
elev = srtm.get_data()

geod = Geod(ellps="WGS84")
all_samples = []
obstacles = []

print("\nProcessing azimuths...\n")

for az in tqdm(range(0, 360, int(az_step))):
    samples = analyze_azimuth(
        lat, lon,
        tower_height,
        freq_mhz,
        az,
        max_dist_km,
        sample_m,
        rx_height,
        fresnel_fraction,
        geod,
        elev
    )
    all_samples.extend(samples)

    for s in samples:
        if s["clearance_m"] < 0 or s["knife_edge_loss_db"] > 0:
            obstacles.append(s)


# Save JSON
with open("output_results.json", "w") as f:
    json.dump(all_samples, f, indent=2)

print("\n✔ Saved calculations as output_results.json")


# -----------------------------
# 3D MAP GENERATION
# -----------------------------
print("\nGenerating 3D terrain model (may take ~1 minute)...")

lats, lons, heights = build_3d_map(lat, lon, max_dist_km, elev)


# Surface
surface = go.Surface(
    z=heights,
    x=lons,
    y=lats,
    showscale=False,
    opacity=0.9
)

# Obstacles
obs_lat  = [o["lat"] for o in obstacles]
obs_lon  = [o["lon"] for o in obstacles]
obs_hgt  = [o["terrain_elev_m"] for o in obstacles]

obstacle_points = go.Scatter3d(
    x=obs_lon,
    y=obs_lat,
    z=obs_hgt,
    mode='markers',
    marker=dict(size=4, color="red"),
    name="Obstacles"
)

# Tower
tower_ground = safe_elev(elev.get_elevation(lat, lon))

tower_marker = go.Scatter3d(
    x=[lon],
    y=[lat],
    z=[tower_ground + tower_height],
    mode='markers',
    marker=dict(size=8, color="blue"),
    name="Tower"
)


layout = go.Layout(
    title="3D Terrain Map with Obstacles",
    scene=dict(
        xaxis_title="Longitude",
        yaxis_title="Latitude",
        zaxis_title="Elevation (m)"
    ),
    margin=dict(l=0, r=0, b=0, t=30)
)

fig = go.Figure(data=[surface, obstacle_points, tower_marker], layout=layout)

output_html = "3d_map.html"
pyo.plot(fig, filename=output_html, auto_open=False)

print(f"\n✔ 3D interactive map saved as: {output_html}")
print("Open this file in any browser (Chrome/Edge).\n")

print("=== DONE ===\n")
