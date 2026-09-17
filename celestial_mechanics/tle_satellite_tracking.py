import numpy as np
import matplotlib.pyplot as plt
from sgp4.api import Satrec, jday
from datetime import datetime
from astropy.coordinates import EarthLocation, AltAz, ICRS
from astropy.time import Time
from astropy import units as u
from math import radians, degrees, atan2, asin, sqrt

# Step 1: Define your location (Tromsø, Norway)
latitude = 69.682793
longitude = 18.9754455
elevation = 0  # Elevation in meters (approximate)

# Convert coordinates to an EarthLocation object in Astropy
location = EarthLocation.from_geodetic(longitude * u.deg, latitude * u.deg, elevation * u.m)


# Step 2: Load the TLE data from the file
def load_tle(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        satellites = []
        for i in range(0, len(lines), 3):
            name = lines[i].strip()  # Satellite name
            line1 = lines[i + 1].strip()  # TLE line 1
            line2 = lines[i + 2].strip()  # TLE line 2
            satellites.append((name, line1, line2))
        return satellites


# Load TLE data from 'TLE.txt'
tle_data = load_tle('TLE.txt')

# Parse the TLE for ASBM-1 and ASBM-2
sat1_name, tle_asbm_1_1, tle_asbm_1_2 = tle_data[0]
sat2_name, tle_asbm_2_1, tle_asbm_2_2 = tle_data[1]

# Step 3: Create satellite objects using SGP4
sat1 = Satrec.twoline2rv(tle_asbm_1_1, tle_asbm_1_2)
sat2 = Satrec.twoline2rv(tle_asbm_2_1, tle_asbm_2_2)

# Step 4: Use the current time for propagation instead of the TLE epoch
current_time = datetime.utcnow()  # Use UTC time for real-time propagation

# Calculate the Julian Date for the current time
jd, fr = jday(current_time.year, current_time.month, current_time.day, current_time.hour, current_time.minute,
              current_time.second)

# Propagate satellites' positions using the current time
e, r, v = sat1.sgp4(jd, fr)
e2, r2, v2 = sat2.sgp4(jd, fr)

# Step 5: Check for propagation errors
if e != 0:
    print(f"Error propagating satellite 1: {e}")
else:
    print(f"Satellite 1 position (Cartesian): {r}")
if e2 != 0:
    print(f"Error propagating satellite 2: {e2}")
else:
    print(f"Satellite 2 position (Cartesian): {r2}")


# Step 6: Convert Cartesian coordinates to latitude/longitude
def cartesian_to_lonlat(x, y, z):
    """Convert Cartesian coordinates to latitude and longitude."""
    r = sqrt(x ** 2 + y ** 2 + z ** 2)  # Distance from Earth's center
    latitude = asin(z / r)  # Latitude in radians
    longitude = atan2(y, x)  # Longitude in radians

    # Convert radians to degrees
    latitude = degrees(latitude)
    longitude = degrees(longitude)

    return latitude, longitude


# Convert satellite positions to latitude and longitude
lat1, lon1 = cartesian_to_lonlat(r[0], r[1], r[2])
lat2, lon2 = cartesian_to_lonlat(r2[0], r2[1], r2[2])


# Adjust the longitude to match conventions (-180 to 180 degrees)
def adjust_longitude(lon):
    """Adjust the longitude for correct hemisphere representation (-180 to 180)."""
    # Ensure longitude is in the -180 to 180 range
    if lon > 180:
        lon -= 360  # Adjust to -180 to 180 range
    elif lon < -180:
        lon += 360  # Adjust to -180 to 180 range
    return lon


# Apply longitude adjustment
lon1 = adjust_longitude(lon1)
lon2 = adjust_longitude(lon2)

# Debug prints to verify longitude and latitude
print(f"Satellite 1: Latitude = {lat1:.2f}, Longitude = {lon1:.2f}")
print(f"Satellite 2: Latitude = {lat2:.2f}, Longitude = {lon2:.2f}")

# Step 7: Convert to ICRS (International Celestial Reference System)
from astropy.coordinates import CartesianRepresentation

satellite1_position = CartesianRepresentation(x=r[0] * u.km, y=r[1] * u.km, z=r[2] * u.km)
satellite2_position = CartesianRepresentation(x=r2[0] * u.km, y=r2[1] * u.km, z=r2[2] * u.km)

# Convert to ICRS (International Celestial Reference System)
icrs1 = ICRS(satellite1_position)
icrs2 = ICRS(satellite2_position)

# Step 8: Convert to AltAz (Altitude/Azimuth) frame based on your location
obs_time = Time(current_time)
altaz_frame = AltAz(obstime=obs_time, location=location)

altaz1 = icrs1.transform_to(altaz_frame)
altaz2 = icrs2.transform_to(altaz_frame)

# Step 9: Plotting the orbit and the satellite locations
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='mollweide')

# Plot the Earth as a circle (no image needed)
earth_circle = plt.Circle((0, 0), 1, color='lightblue', ec='black', lw=1)
ax.add_artist(earth_circle)

# Plot the satellite locations on the Mollweide projection
ax.scatter(radians(lon1), radians(lat1), c='red', marker='x', label=sat1_name)
ax.scatter(radians(lon2), radians(lat2), c='blue', marker='x', label=sat2_name)

# Plot your location on Earth (Tromsø)
ax.scatter(radians(longitude), radians(latitude), color="green", marker="o", label="Your Location")

# Show altitude and azimuth for both satellites
print(f"{sat1_name} Altitude: {altaz1.alt:.2f}, Azimuth: {altaz1.az:.2f}")
print(f"{sat2_name} Altitude: {altaz2.alt:.2f}, Azimuth: {altaz2.az:.2f}")

# Show the plot
ax.legend()
plt.title('Satellite Positions and Your Location')
plt.show()

