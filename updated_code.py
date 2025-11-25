import math
import json
import numpy as np
from pyproj import Geod
import srtm
from tqdm import tqdm
import plotly.graph_objs as go
import plotly.offline as pyo

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
    # Classic approximate formula
    if v <= -0.78:
        return 0.0
    return 6.9 + 20 * math.log10(math.sqrt((v - 0.1) ** 2 + 1) + v - 0.1)

def safe_elev(e):
    return e if e is not None else 0.0

# -----------------------------
# Fresnel radius (using freq in Hz)
# -----------------------------
def fresnel_radius_freq(freq_hz, d1, d2):
    if d1 + d2 == 0:
        return 0.0
    lam = 299792458.0 / freq_hz
    return math.sqrt((lam * d1 * d2) / (d1 + d2))

# -----------------------------
# Earth bulge for a point at distance d along path of length D
# -----------------------------
def earth_bulge(d, D):
    # d*(D-d)/(2*R_eff)
    return d * (D - d) / (2.0 * R_EFF)

# -----------------------------
# Core coverage computation using bulge + Fresnel approach
# -----------------------------
def compute_coverage(lat, lon, height_m, freq_mhz,
                     az_step_deg=1, max_dist_km=10, sample_m=30,
                     rx_h=1.5, fresnel_frac=0.6, sensitivity_db=25):
    """
    Returns per-azimuth max usable distance (LOS or diffraction <= sensitivity_db).
    - az_step_deg: 1 for every degree
    - sample_m: sampling resolution (meters)
    - fresnel_frac: fraction of first Fresnel zone to require clear (e.g. 0.6)
    - sensitivity_db: loss threshold above which link is considered unusable
    """
    freq_hz = freq_mhz * 1e6
    geod = Geod(ellps="WGS84")

    max_dist = max_dist_km * 1000.0
    # Start from sample_m (avoid distance 0)
    dists = np.arange(sample_m, max_dist + sample_m, sample_m)

    dem = srtm.get_data()
    tower_ground = safe_elev(dem.get_elevation(lat, lon))
    tower_abs = tower_ground + height_m

    results = []
    full_samples = []

    for az in tqdm(range(0, 360, int(az_step_deg))):
        # We'll find the farthest D (from dists) such that no intermediate obstruction
        # causes loss > sensitivity_db when evaluating the path to that D.
        max_usable_D = 0.0
        max_usable_mode = "LOS"
        max_usable_loss = 0.0
        max_usable_coord = (lat, lon)

        blocked = False

        # iterate candidate receiver distances D (increasing)
        for i, D in enumerate(dists):
            # receiver coord at distance D
            lonR, latR, _ = geod.fwd(lon, lat, az, D)
            elevR = safe_elev(dem.get_elevation(latR, lonR))
            rx_abs = elevR + rx_h

            # keep track whether this D is usable
            D_usable = True
            D_max_loss = 0.0
            D_mode = "LOS"

            # check all intermediate points along path (from near to D)
            # we exclude the receiver point itself from intermediate checks
            for j in range(i):
                d = dists[j]
                lonI, latI, _ = geod.fwd(lon, lat, az, d)
                elevI = safe_elev(dem.get_elevation(latI, lonI))

                # line height at intermediate point (straight line between tx_abs and rx_abs)
                line_h = tower_abs + (rx_abs - tower_abs) * (d / D)

                # include earth bulge
                bulge = earth_bulge(d, D)

                # distances for Fresnel calc
                d1 = d
                d2 = D - d
                if d2 <= 0 or (d1 + d2) == 0:
                    continue

                # first Fresnel radius (meters)
                f1 = fresnel_radius_freq(freq_hz, d1, d2)
                clearance_required = f1 * fresnel_frac

                # obstacle top (terrain + bulge)
                obst = elevI + bulge

                # if obstacle intrudes into the required clearance zone
                if obst >= (line_h - clearance_required):
                    # penetration height
                    h = obst - line_h
                    # compute v parameter for knife-edge (use lam)
                    lam = 299792458.0 / freq_hz
                    # guard denom
                    denom = lam * d1 * d2
                    if denom == 0:
                        v = 0.0
                    else:
                        v = h * math.sqrt(2.0 * (d1 + d2) / denom)

                    loss = knife_edge_loss_db(v)
                    D_max_loss = max(D_max_loss, loss)

                    # if loss exceeds sensitivity, this D is unusable
                    if loss > sensitivity_db:
                        D_usable = False
                        D_mode = "Blocked"
                        break
                    else:
                        # still usable (diffraction) but mark mode
                        D_mode = "Diffraction"
                        # continue checking other intermediate points for worse loss

            # store detailed sample for output if needed
            full_samples.append({
                "azimuth": az,
                "distance_m": D,
                "lat": latR,
                "lon": lonR,
                "terrain_elev_m": elevR,
                "rx_abs_m": rx_abs,
                "usable": D_usable,
                "mode": D_mode,
                "max_loss_db_on_path": D_max_loss
            })

            if not D_usable:
                # first distance that fails → stop here, max usable is previous distance
                # if i==0 then nothing was usable (max_usable_D remains 0)
                blocked = True
                break
            else:
                # this D is usable (LOS or acceptable diffraction) → record as current max
                max_usable_D = D
                max_usable_mode = D_mode
                max_usable_loss = D_max_loss
                max_usable_coord = (latR, lonR)

        # if we never found blocking up to max_dist, treat as LOS to max_dist
        if not blocked and max_usable_D >= max_dist - 1e-6:
            max_usable_mode = "LOS"
            max_usable_loss = 0.0

        results.append({
            "azimuth": az,
            "max_distance_m": max_usable_D,
            "max_distance_km": max_usable_D / 1000.0,
            "coord_lat": max_usable_coord[0],
            "coord_lon": max_usable_coord[1],
            "mode": max_usable_mode,
            "loss_db": max_usable_loss
        })

    return results, full_samples

# -----------------------------
# 3D Terrain builder (unchanged)
# -----------------------------
def build_3d_map(center_lat, center_lon, radius_km, elev):
    GRID = 200
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
# MAIN (user I/O)
# =================================================================
if __name__ == "__main__":
    print("\n=== Tower Obstacle Detection + Directional Max Coverage ===\n")

    lat = float(input("Enter tower latitude  : "))
    lon = float(input("Enter tower longitude : "))
    tower_height = float(input("Enter tower height (m) : "))
    freq_mhz = float(input("Enter frequency (MHz)  : "))
    az_step = int(float(input("Azimuth step (degrees, default 1): ") or 1))
    max_dist_km = float(input("Max distance km (default 10 km): ") or 10.0)
    sample_m = float(input("Sampling distance meters (default 30m): ") or 30.0)
    rx_height = float(input("Receiver height (default 1.5m): ") or 1.5)
    fresnel_fraction = float(input("Fresnel fraction (default 0.6): ") or 0.6)
    sensitivity_db = float(input("Max acceptable knife-edge loss (dB) (default 25): ") or 25.0)

    print("\nLoading SRTM elevation data...")
    elev = srtm.get_data()

    print("\nComputing coverage per azimuth (this can take a while)...\n")
    results, full_samples = compute_coverage(
        lat, lon, tower_height, freq_mhz,
        az_step_deg=az_step, max_dist_km=max_dist_km, sample_m=sample_m,
        rx_h=rx_height, fresnel_frac=fresnel_fraction, sensitivity_db=sensitivity_db
    )

    # print & save per-degree summary
    with open("coverage_direction_report.json", "w") as f:
        json.dump(results, f, indent=2)

    with open("output_results.json", "w") as f:
        json.dump(full_samples, f, indent=2)

    print("\n=== Directional Coverage Summary ===\n")
    for r in results:
        lat_r = r["coord_lat"]
        lon_r = r["coord_lon"]
        km = r["max_distance_km"]
        print(f"Azimuth {r['azimuth']:3d}°:  {km:6.3f} km  → ({lat_r:.6f}, {lon_r:.6f})  Mode: {r['mode']}  Loss: {r['loss_db']:.2f} dB")

    print("\n✔ Saved direction coverage as coverage_direction_report.json")
    print("✔ Saved full samples as output_results.json")

    # -----------------------------
    # 3D map generation (optional)
    # -----------------------------
    create_map = input("\nGenerate 3D terrain map (y/N)? ").lower().startswith("y")
    if create_map:
        print("Generating 3D terrain model (this may take ~1 minute)...")
        lats, lons, heights = build_3d_map(lat, lon, max_dist_km, elev)

        surface = go.Surface(
            z=heights,
            x=lons,
            y=lats,
            showscale=False,
            opacity=0.9
        )

        # mark obstacles from full_samples where usable==False (points that caused blocking)
        blocked_points = [s for s in full_samples if s.get("usable") is False]
        obs_lat = [s["lat"] for s in blocked_points]
        obs_lon = [s["lon"] for s in blocked_points]
        obs_hgt = [s["terrain_elev_m"] for s in blocked_points]

        obstacle_points = go.Scatter3d(
            x=obs_lon,
            y=obs_lat,
            z=obs_hgt,
            mode='markers',
            marker=dict(size=4, color="red"),
            name="Blocking Points"
        )

        tower_ground = safe_elev(elev.get_elevation(lat, lon))
        tower_marker = go.Scatter3d(
            x=[lon],
            y=[lat],
            z=[tower_ground + tower_height],
            mode='markers',
            marker=dict(size=8, color="blue"),
            name="Tower"
        )

        fig = go.Figure(data=[surface, obstacle_points, tower_marker])
        fig.update_layout(
            title="3D Terrain Map with Blocking Points",
            scene=dict(xaxis_title="Longitude", yaxis_title="Latitude", zaxis_title="Elevation (m)"),
            margin=dict(l=0, r=0, b=0, t=30)
        )

        output_html = "3d_map.html"
        pyo.plot(fig, filename=output_html, auto_open=False)
        print(f"\n✔ 3D interactive map saved as: {output_html}")

    print("\n=== DONE ===\n")
