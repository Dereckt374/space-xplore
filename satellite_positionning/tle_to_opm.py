from skyfield.api import EarthSatellite, load, wgs84, N, S, E, W
import numpy as np
from math import degrees

# from astropy import units as u

# from poliastro.bodies import Earth, Mars, Sun
# from poliastro.twobody import Orbit


def TLE_to_orbital_parameters(TLE):
    """
    Converti une TLE en paramètres orbitaux
    formules :
        $$ \frac{a^3}{P^2} =\frac{\mu}{4 \pi^2} $$
        $$  a^3 =\frac{\mu}{4 \pi^2} \times P^2 $$    
        ![](https://wikimedia.org/api/rest_v1/media/math/render/svg/8ccc0844c5cdd7b0a9375206336284d919aa33a8)
    """

    mu = 3.9860044188e14 #[m3/s2]
    Rt = 6371 #km
    list_param = {}
    L1, L2 = TLE.splitlines()
    sat = EarthSatellite(L1, L2)
    mean_mo = sat.model.no_kozai # no_kozai: mean motion (radians/minute)
    list_param['mean_motion_deg_per_min'] = degrees(mean_mo)
    P_orbit = 2 * np.pi / mean_mo * 60   #min
    list_param['orbital_period_min'] = P_orbit / 60
    h = (mu / (4 * np.pi**2) * P_orbit ** 2 )**(1/3) / 1000 - Rt   ## Altitude of an orbit, km
    list_param['orbit_altitude_km'] = h
    Ma = sat.model.mo # mo: mean anomaly (radians)
    e = sat.model.ecco  # ecco: eccentricity
    list_param['mean_anomaly_deg'] = degrees(Ma)
    list_param['eccentricity'] = e
    tru_anomaly = Ma + (2*e - 1/4 *e**3) * np.sin(Ma) + 5/4 * e**2 * np.sin(2*Ma) + 13/12 * e**3 * np.sin(3*Ma) #rad
    list_param['tru_anomaly_deg'] = degrees(tru_anomaly)
    inc = sat.model.inclo  # inclo: inclination (radians)
    raan = sat.model.nodeo        # nodeo: right ascension of ascending node (radians)
    arg_per = sat.model.argpo        # argpo: argument of perigee (radians)
    list_param['inclination_deg'] = degrees(inc)
    list_param['raan_deg'] = degrees(raan)
    list_param['arg_perigee_deg'] = degrees(arg_per)
    list_param['B_star'] = sat.model.bstar
    return(list_param)