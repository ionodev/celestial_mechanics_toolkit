import numpy as np

def drag_on_sattelite(altitude):
    # For this exercise I've assumed a 3U cubesat for the spacecraft
    m = 4 #kg spacecraft mass
    S = 0.1562 #m^2 spacecraft surface area'
    G = 6.6743e-11 # m^3 / (kg s^2) Gravitational constant
    M = 5.972e24 #kg Earth mass
    R = altitude*1000 + 6371e3 # orbit radius including earth radius and altitude
    T = np.sqrt(4*np.pi**2 * R**3 / (G*M)) # s period (kepler 3 law)
    C_D = 2.5 # drag coefficient (usually between 2 and 4)

    def mass_density(target_altitude, file_path='MSIS.dat'):
        """
        Given a target altitude in km, this function reads the MSIS data from a file,
        and returns the mass density at the specified altitude in kg/m^3. The altitude must be in the data file.
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
        densities = np.array(densities)

        # Check if the target altitude exists in the dataset
        if target_altitude in heights:
            # Return the mass density for the exact altitude, converted to kg/m^3
            index = np.where(heights == target_altitude)[0][0]
            mass_density_kg_per_m3 = densities[index] * 1000  # Convert g/cm^3 to kg/m^3
            return mass_density_kg_per_m3
        else:
            raise ValueError(f"Altitude {target_altitude} km is not available in the data.")

    d_T = T * -3*np.pi*mass_density(altitude)*R*C_D*(S/m)
    orbit_frac_percent = np.abs(d_T)/T*100
    print(f'@ Altitude: {altitude} km | Orbital period T = {T:.3f} s | ΔT = {d_T:.3f} s/day | ΔT/T = {orbit_frac_percent:.3f} %')

drag_on_sattelite(300) #km
drag_on_sattelite(200) #km
drag_on_sattelite(150) #km
drag_on_sattelite(125) #km
drag_on_sattelite(120) #km
drag_on_sattelite(115) #km
drag_on_sattelite(110) #km
drag_on_sattelite(105) #km
drag_on_sattelite(100) #km
