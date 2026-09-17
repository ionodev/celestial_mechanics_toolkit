import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# constants
G = 6.6743e-11  # Gravitational constant
A = 0.1562  # m^2 (cross-sectional area)
m = 4  # kg (mass)
M = 5.97219e24  # kg (mass of Earth)
R_e = 6371  # km (Earth radius)
C_d = 2.5  # Drag coefficient
d_t = 60  # 1 minute time step (start with a small value)


# Exponential model for density vs. altitude
def exponential_decay(h, rho_0, H):
    return rho_0 * np.exp(-h / H)


def mass_density(target_altitude, file_path='MSIS.dat'):
    """
    Given a target altitude in km, this function reads the MSIS data from a file,
    and returns the mass density at the specified altitude in kg/m^3.
    If the altitude is above 600 km, the function extrapolates using the exponential decay model.
    """
    heights = []
    densities = []

    # Read the MSIS data file
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('%'):
                continue  # Skip comment lines

            values = line.split()
            if len(values) < 6:
                continue  # Skip malformed lines

            height = float(values[0])  # Height in km
            mass_density_g_per_cm3 = float(values[4])  # Mass density in g/cm^3

            heights.append(height)
            densities.append(mass_density_g_per_cm3)

    heights = np.array(heights)
    densities = np.array(densities) * 1000  # Convert g/cm^3 to kg/m^3

    # If the target altitude is greater than 600 km, extrapolate using the exponential decay model
    if target_altitude <= 600:
        # Use np.interp() to interpolate the mass density for the target altitude
        interpolated_density = np.interp(target_altitude, heights, densities)
    else:
        # Fit the exponential decay model to the data (up to 600 km)
        params, _ = curve_fit(exponential_decay, heights, densities, p0=[1e-3, 700])  # Adjust the initial guess
        rho_0, H = params

        # Extrapolate, but reduce the decay rate above 600 km
        interpolated_density = exponential_decay(target_altitude, rho_0, H)

        # To prevent densities from becoming too small, set a lower threshold for density
        # We choose a reasonable floor value for the density at high altitudes
        if interpolated_density < 1e-10:
            interpolated_density = 1e-10  # Minimum threshold for density at 1000 km

    return interpolated_density


def alpha_d(p, V):
    return -0.5 * p * C_d * (A / m) * V ** 2


def V(R):
    return np.sqrt(G * M / R)


def T(R):
    return 2 * np.pi * np.sqrt(R ** 3 / (G * M))


def simulate_trajectory(initial_altitude):
    height = initial_altitude
    t = 0  # initial time
    time = [t]
    heights = [height]

    while height >= 200:  # until height is 200km
        R = R_e + height
        p = mass_density(height)  # Get mass density
        v = V(R)  # Velocity
        drag_force = alpha_d(p, v)  # Drag force

        # Debug: Print relevant variables to check the issue
        if t % (60 * 60) == 0:  # Print every hour
            print(
                f"Time: {t / 3600:.2f} hours, Height: {height:.2f} km, Density: {p:.2e}, Velocity: {v:.2f} km/s, Drag Force: {drag_force:.2e}")

        d_R = drag_force * T(R) / np.pi * d_t  # Change in height due to drag
        height += d_R
        t += d_t
        heights.append(height)
        time.append(t)

    return np.array(time), np.array(heights)  # Return time in seconds, heights in km


# Simulate trajectories for an initial altitude of 1000 km
time_1000, heights_1000 = simulate_trajectory(1000)

# Plot: 1000 km initial altitude (Time in days)
fig, ax = plt.subplots()

ax.plot(time_1000 / (3600 * 24), heights_1000, label='Initial Altitude = 1000 km', color='g')  # Time in days
ax.set_xlabel("Time (days)")
ax.set_ylabel("Height (km)")
ax.set_title("Spacecraft Altitude vs. Time (Initial Altitude = 1000 km)")
ax.legend(loc="upper right")
ax.grid(True)

plt.show()
