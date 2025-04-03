import math
from datetime import datetime
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), '..')))
import fct_space 

# Constante gravitationnelle de la Terre (km^3/s^2)
MU = 398600.4418


# Exemple d'utilisation
if __name__ == "__main__":
    raan = 120.0  # RAAN en degrés
    inclination = 97.5  # Inclinaison en degrés
    arg_perigee = 45.0  # Argument du périgée en degrés
    sma = 6871.0  # Semi-major axis en km pour une orbite basse (~500 km altitude)
    eccentricity = 0.0012  # Excentricité (typique pour cubesat en LEO)
    date_time = datetime(2024, 10, 2, 12, 0, 0)  # Date/Time (UTC)

    # Générer TLE
    tle_line_1, tle_line_2 = fct_space.generate_tle(
        raan, inclination, arg_perigee, sma, eccentricity, date_time
    )

    print("TLE Generated:")
    print(tle_line_1)
    print(tle_line_2)
