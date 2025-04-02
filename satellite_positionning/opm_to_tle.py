import math
from datetime import datetime

# Constante gravitationnelle de la Terre (km^3/s^2)
MU = 398600.4418


def mean_motion_from_sma(sma_km):
    """Calculer la Mean Motion (révolution par jour) à partir du Semi-Major Axis (km)"""
    # Mean Motion en rad/s
    mean_motion_rad_s = math.sqrt(MU / (sma_km**3))

    # Conversion en révolutions par jour (1 tour = 2 * pi rad)
    mean_motion_rev_day = mean_motion_rad_s * 86400 / (2 * math.pi)

    return mean_motion_rev_day


def format_tle_element(element, width, decimals=0):
    """Format un élément TLE avec largeur et décimales"""
    return f"{element:0{width}.{decimals}f}"


def create_tle_line_1(
    satellite_number,
    classification,
    launch_year,
    launch_number,
    launch_piece,
    epoch_year,
    epoch_day,
):
    """Créer la première ligne TLE"""
    line_1 = f"1 {satellite_number:05}{classification} {launch_year:02}{launch_number:03}{launch_piece} {epoch_year:02}{format_tle_element(epoch_day, 12, 8)}  .00000000  00000-0  00000-0 0  0000"
    return line_1


def create_tle_line_2(
    satellite_number,
    inclination,
    raan,
    eccentricity,
    arg_perigee,
    mean_anomaly,
    mean_motion,
):
    """Créer la deuxième ligne TLE"""
    eccentricity_formatted = f"{int(eccentricity * 1e7):07d}"
    line_2 = f"2 {satellite_number:05} {format_tle_element(inclination, 8, 4)} {format_tle_element(raan, 8, 4)} {eccentricity_formatted} {format_tle_element(arg_perigee, 8, 4)} {format_tle_element(mean_anomaly, 8, 4)} {format_tle_element(mean_motion, 11, 8)}00000"
    return line_2


def date_to_epoch_days(dt):
    """Convertir une date en jour fractionnaire pour TLE"""
    year = dt.year
    # Si année après 2000, ajuster pour le TLE
    epoch_year = year % 100
    # Jour de l'année + fraction du jour
    epoch_day = (
        dt.timetuple().tm_yday + dt.hour / 24 + dt.minute / 1440 + dt.second / 86400
    )
    return epoch_year, epoch_day


def generate_tle(raan, inclination, arg_perigee, sma, eccentricity, dt):
    """Générer un TLE approximatif à partir de paramètres orbitaux"""

    # Satellite et lancement fictif (cubesat 3U)
    satellite_number = 99999
    classification = "U"  # Unclassified
    launch_year = 23  # 2023
    launch_number = 99  # Lancement fictif
    launch_piece = "A"  # Fragment A

    # Calculer la Mean Motion (révolutions par jour) à partir du SMA
    mean_motion = mean_motion_from_sma(sma)

    # Date et conversion en jour TLE (epoch)
    epoch_year, epoch_day = date_to_epoch_days(dt)

    # Anomalie moyenne fictive (peut être dérivée avec plus de données)
    mean_anomaly = 0.0  # Supposons 0 pour simplifier

    # Créer les lignes TLE
    line_1 = create_tle_line_1(
        satellite_number,
        classification,
        launch_year,
        launch_number,
        launch_piece,
        epoch_year,
        epoch_day,
    )
    line_2 = create_tle_line_2(
        satellite_number,
        inclination,
        raan,
        eccentricity,
        arg_perigee,
        mean_anomaly,
        mean_motion,
    )

    return line_1, line_2


# Exemple d'utilisation
if __name__ == "__main__":
    raan = 120.0  # RAAN en degrés
    inclination = 97.5  # Inclinaison en degrés
    arg_perigee = 45.0  # Argument du périgée en degrés
    sma = 6871.0  # Semi-major axis en km pour une orbite basse (~500 km altitude)
    eccentricity = 0.0012  # Excentricité (typique pour cubesat en LEO)
    date_time = datetime(2024, 10, 2, 12, 0, 0)  # Date/Time (UTC)

    # Générer TLE
    tle_line_1, tle_line_2 = generate_tle(
        raan, inclination, arg_perigee, sma, eccentricity, date_time
    )

    print("TLE Generated:")
    print(tle_line_1)
    print(tle_line_2)
