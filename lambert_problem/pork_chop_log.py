import numpy as np
import matplotlib.pyplot as plt
import pykep as pk
from tqdm import tqdm
import numpy.ma as ma

# === Constants ===
MU_SUN = pk.MU_SUN
t_impact = 2200  # days after MJD2000

# === Define Earth ===
earth = pk.planet.jpl_lp("earth")

# === Define impact epoch ===
t_impact_epoch = pk.epoch(t_impact, 'mjd2000')

# === Get Earth's position at impact (where the comet will hit) ===
r_impact, v_impact = earth.eph(t_impact_epoch)

# === Comet's known velocity vector (given) ===
v_comet = np.array([-15e3, -28e3, 28e3])

# === Define the comet's orbit ===
comet = pk.planet.keplerian(
    t_impact_epoch,   # epoch
    r_impact,         # position vector (m)
    v_comet,          # velocity vector (m/s)
    MU_SUN,           # gravitational parameter of Sun
    1e3,              # radius (m)
    1e3,              # safe radius (m)
    1e12,             # mass (kg)
    "Comet"           # name
)

# === Verify comet orbit ===
print("\n🪐 Comet Keplerian Elements:")
a, e, i, RAAN, arg_peri, M = comet.osculating_elements(t_impact_epoch)
print(f"  a = {a/1.496e11:.3f} AU, e = {e:.4f}, i = {np.degrees(i):.2f}°")

# === Pork-Chop Parameters ===
dep_start = 0
dep_end = 2250
arr_start = 100
arr_end = 2300 - 5
n_points = 2000  # grid resolution

departure_dates = np.linspace(dep_start, dep_end, n_points)
arrival_dates = np.linspace(arr_start, arr_end, n_points)
DEP, ARR = np.meshgrid(departure_dates, arrival_dates)
DELTA_V = np.full_like(DEP, np.nan, dtype=float)

# === Compute Pork-Chop Grid ===
print("\n🛠 Computing pork-chop grid...")
for i in tqdm(range(n_points), desc="Rows"):
    for j in range(n_points):
        t_dep = DEP[i, j]
        t_arr = ARR[i, j]
        if t_arr <= t_dep:
            continue  # invalid transfer
        epoch_dep = pk.epoch(t_dep, 'mjd2000')
        epoch_arr = pk.epoch(t_arr, 'mjd2000')
        try:
            r1, v1 = earth.eph(epoch_dep)
            r2, v2 = comet.eph(epoch_arr)
            dt = (t_arr - t_dep) * pk.DAY2SEC
            lambert_sol = pk.lambert_problem(r1, r2, dt, MU_SUN)
            v_dep = np.array(lambert_sol.get_v1()[0])
            delta_v = np.linalg.norm(v_dep - np.array(v1))
            if 0 < delta_v < 1e6:  # filter nonsense
                DELTA_V[i, j] = delta_v
        except Exception:
            continue

# === Ensure DELTA_V is populated ===
if np.all(np.isnan(DELTA_V)):
    print("Error: DELTA_V was not populated. Check the pork-chop grid calculation.")
else:
    print(f"✅ Pork-chop grid computed. Valid Δv range: {np.nanmin(DELTA_V):.1f} – {np.nanmax(DELTA_V):.1f} m/s")

# === Mask Δv outside 0-90 km/s ===
dv_min = 0
dv_max = 50e3  # 90 km/s in m/s
masked_dv_linear = ma.masked_invalid(DELTA_V)
masked_dv_linear = ma.masked_outside(masked_dv_linear, dv_min, dv_max)

# === Convert Δv to km/s for plotting ===
masked_dv_kms = masked_dv_linear / 1000.0  # m/s -> km/s
dv_min_kms = dv_min / 1000.0
dv_max_kms = dv_max / 1000.0

# === Plotting ===
plt.figure(figsize=(13, 6))
img = plt.pcolormesh(
    DEP,
    ARR,
    masked_dv_kms,
    vmin=dv_min_kms,
    vmax=dv_max_kms,
    shading='auto',
    cmap='viridis'
)

cbar = plt.colorbar(img)
cbar.set_label('Δv (km/s)')  # Colorbar in km/s

plt.xlabel("Departure Date (days after MJD2000)")
plt.ylabel("Arrival Date (days after MJD2000)")
plt.title("🚀 Pork-Chop Plot: Earth-to-Comet Transfer Δv (linear scale, km/s)")

# === Vertical lines at 5, 10, 20, 40, 365 days before impact ===
critical_times = [5, 10, 20, 40, 365]
for delta in critical_times:
    dep_time = t_impact - delta
    plt.axvline(x=dep_time, color='red', linestyle='--', alpha=0.7)
    plt.text(dep_time + 2, 2110, f'-{delta} d', color='red', rotation=90)

# === Horizontal impact line ===
plt.axhline(y=t_impact, color='white', linestyle='-', linewidth=2, alpha=0.9)
plt.text(1810, t_impact + 5, "Comet Impact (t = 2200)", color='white')

# === Print Δv values at the intersections to console in m/s (original values) ===
for delta in critical_times:
    dep_time = t_impact - delta
    arr_time = t_impact

    i = np.argmin(np.abs(arrival_dates - arr_time))
    j = np.argmin(np.abs(departure_dates - dep_time))

    dv_value = masked_dv_linear[i, j]  # m/s
    print(f"Δv at departure = {dep_time:.1f}, arrival = {arr_time:.1f}  -->  {dv_value:.1f} m/s")

# === Set axis limits ===
plt.xlim(1800, 2250)
plt.ylim(2100, 2300)

plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
