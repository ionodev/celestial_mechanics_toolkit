import numpy as np
import matplotlib.pyplot as plt
import pykep as pk
from tqdm import tqdm
import numpy.ma as ma
from scipy.ndimage import gaussian_filter

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
dep_end = 2100
arr_start = 100
arr_end = 2300 - 5
n_points = 200  # Higher resolution grid

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
            continue
        epoch_dep = pk.epoch(t_dep, 'mjd2000')
        epoch_arr = pk.epoch(t_arr, 'mjd2000')
        try:
            r1, v1 = earth.eph(epoch_dep)
            r2, v2 = comet.eph(epoch_arr)
            dt = (t_arr - t_dep) * pk.DAY2SEC
            lambert_sol = pk.lambert_problem(r1, r2, dt, MU_SUN)
            v_dep = np.array(lambert_sol.get_v1()[0])
            delta_v = np.linalg.norm(v_dep - np.array(v1))
            if 0 < delta_v < 1e6:
                DELTA_V[i, j] = delta_v
        except Exception:
            continue

# === Mask invalid values ===
masked_dv = ma.masked_invalid(DELTA_V)
masked_dv = ma.masked_less_equal(masked_dv, 0)

# === Optional smoothing for better visuals ===
smooth_dv = gaussian_filter(masked_dv, sigma=1.0)

# === Plot (linear scale, limited range) ===
vmin, vmax = 6e3, 1.5e4
plt.figure(figsize=(14, 7), dpi=200)
img = plt.pcolormesh(
    DEP,
    ARR,
    smooth_dv,
    vmin=vmin,
    vmax=vmax,
    shading='auto',
    cmap='viridis'
)

cbar = plt.colorbar(img)
cbar.set_label('Δv (m/s)')
cbar.ax.tick_params(which='both', direction='in')

plt.xlabel("Departure Date (days after MJD2000)")
plt.ylabel("Arrival Date (days after MJD2000)")
plt.title("🚀 Pork-Chop Plot: Earth-to-Comet Transfer Δv (linear scale)")

# === Add reference lines ===
critical_times = [5, 10, 20, 40, 365]
for delta in critical_times:
    plt.axhline(y=t_impact - delta, color='red', linestyle='--', alpha=0.7)
    plt.text(dep_end - 200, t_impact - delta + 5, f'-{delta}d', color='red')

plt.ylim(1800, 2300)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
