from math import cos, sin, sqrt, radians, log10
from skyfield.api import EarthSatellite, load, wgs84, N, S, E, W, Topos
import numpy as np
import math 
from math import degrees
import requests
from datetime import datetime, timedelta, date
from pytz import timezone, all_timezones
import pandas as pd

# from poliastro.bodies import Earth, Mars, Sun
# from poliastro.twobody import Orbit


G=6.67e-11 #cst grav. [m3/kg/s2]
Mt=5.97e24 #masse terre [kg]
mu=G*Mt #[m3/s2]
Rt = 6371000 #Rayon de la terre [m]

# Coords MTP
lat_mtp = 43.637
long_mtp = 3.843

#Paramètre ellipse

hp=200000
ha=35800000
a=(hp+ha+2*Rt)/2
e=0 #exentricité
# a=200 #demigrand axe
# b=a*(1-e**2)
# c=e*a #demi-distance focale
# S=math.pi*a*b



def vitCircul(nu,r):
    return round((nu/r)**(1/2),3)

def periodOrbit(nu, a):
    return round(2*math.pi*((a**3)/nu)**(1/2),3)

def demiGdaxe(nu,T):
    return round((nu*T*T/(4*math.pi**2))**(1/3),3)

def compare_sats(sat1, sat2):
    for attr in sorted(dir(sat1.model)):
        try:
            value1 = getattr(sat1.model, attr)
        except AttributeError:
            continue
        if not isinstance(value1, (str, int, float)):
            continue
        value2 = getattr(sat2.model, attr)
        verdict = '==' if (value1 == value2) else '!='
        print('{:18} {!r:24} {} {!r:24}'.format(attr, value1, verdict, value2))

# def sat_to_orb(sat):
#     mu = 3.9860044188e14 #[m3/s2]
#     Rt = 6371 #km

#     mean_mo = sat.model.no_kozai
#     P_orbit = 2 * np.pi / mean_mo * 60 #min

#     h = (mu / (4 * np.pi**2) * P_orbit ** 2 )**(1/3) / 1000 - Rt   ## Altitude of an orbit
#     Ma = sat.model.mo
#     e = sat.model.ecco
#     tru_anomaly = Ma + (2*e - 1/4 *e**3) * np.sin(Ma) + 5/4 * e**2 * np.sin(2*Ma) + 13/12 * e**3 * np.sin(3*Ma)
#     r_earth = Earth.R_mean << u.km
#     h_orbit = h << u.km 
#     a =  r_earth + h_orbit
#     ecc = e << u.one  ### not used ?
#     inc = degrees(sat.model.inclo) << u.deg
#     raan = degrees(sat.model.nodeo) << u.deg
#     argp = degrees(sat.model.argpo) << u.deg
#     nu = degrees(tru_anomaly) << u.deg
#     orb = Orbit.from_classical(Earth, a, ecc, inc, raan, argp, nu, sat.epoch.to_astropy())
#     return(orb)

def slant_range(theta:float):
    """
    input : elevation, degree
    return le slant range (distance en "vol-d'oiseau") du satellite pour une certaine elevation 'theta'    
    """
    Re = 6378.136 # km, rayon de la terre
    r_orbite_mean = 6928.14 #km

    return ( Re * (    sqrt(r_orbite_mean**2 / Re**2  - cos(radians(theta)) **2  ) - sin(radians(theta))   )     )

def path_loss(S, f):
    """
    return l'atténuation en db de l'atm
    """
    l = 1/f
    return(22 + 20*log10 (S/l))


def TLE_to_sat(tle, name=None):
    """
    Convert a TLE to a EarthSatellite object
    """
    L1, L2 = tle.splitlines()
    sat = EarthSatellite(L1, L2, name)
    return(sat)



def get_tle_from_satnogs(norad_cat_id = 58470):
    r = requests.get('https://db.satnogs.org/api/tle/')       #, auth=(DDP_USERNAME, DDP_PASSWORD))
    try : 
        data = r.json()
        df = pd.DataFrame.from_dict(data)
    except : pass
    df_tle_sat = df[df['norad_cat_id'] == norad_cat_id]
    tle = []
    for i in range(1,3):
        tle.append(df_tle_sat.iloc[:,i].values[0])
    TLE_satnogs = '\n'.join(tle)
    print(f'TLE :\n{TLE_satnogs}')
    return(TLE_satnogs)


def get_all_pass_coord_with_t_beacon(t, TLE):
    """
    t : heure de reception du beacon
    TLE : TLE du satellite à étudier
    output : az et elev en degree
    retourne une trajection du satellite à une certaine date
    """
    satellite = TLE_to_sat(TLE)
    ts = load.timescale()
    t = ts.utc(t.year, t.month, t.day, t.hour,t.minute , range(-10*60,20*60)) # heure utc
    csum = wgs84.latlon(43.637 * N, 3.843 * E) 
    pos = (satellite - csum).at(t)
    alt, az, distance = pos.altaz()
    elev = [round((90 - i), 2) for i in alt.degrees]
    # elev = [round((i), 2) for i in alt.degrees]
    az = [round(i, 2) for i in az.degrees]
    return(az, elev, t)

def polar_coord_beacon(TLE, t_dt):
    """
    TLE : TLE à utiliser pour le passage
    t_dt : liste de datetime correspondant à la reception de chaque beacon
    retourne les coordonnées polaires de reception de chaque beacon
    """
    ts = load.timescale()
    t = ts.from_datetimes([timezone('UTC').localize(i) for i in t_dt])
    satellite = TLE_to_sat(TLE)
    csum = wgs84.latlon(43.637 * N, 3.843 * E)
    pos = (satellite - csum).at(t)
    alt, az, distance = pos.altaz()
    elev = [90 - i for i in alt.degrees]
    return(az.degrees, elev)


def get_last_pass(df, name_col_time='timestamp', start_pass_index = 0,time_max_between_beacon = 20):
    """
    Permet de cropper un df sur seulement un passage, à un index spécifié
    """
    df_last_pass = df[df.index==start_pass_index]
    t_last_received = df[name_col_time][start_pass_index]

    for index, row in df.iterrows():
        if ((t_last_received - df[name_col_time][index]) < pd.Timedelta(time_max_between_beacon, 'minutes')) & ((t_last_received - df[name_col_time][index]) > pd.Timedelta(0, 'seconds')) :
            t_last_received = df[name_col_time][index]
            df_last_pass = pd.concat([df_last_pass, df.iloc[index].to_frame().T])

    number_beacon_decoded = df_last_pass.shape[0]
    df_last_pass = df_last_pass.reset_index(drop=True).convert_dtypes()
    print(f"plage de donnés : {df_last_pass[name_col_time].min().strftime('%d-%m-%Y %H:%M:%S') } - {df_last_pass[name_col_time].max().strftime('%d-%m-%Y %H:%M:%S') }  : {(df_last_pass[name_col_time].max() - df_last_pass[name_col_time].min())}")
    print(f'nombre de beacon lors du passage du {t_last_received.strftime(" %d/%m à %Hh%M")} : {number_beacon_decoded}')
    return(df_last_pass, number_beacon_decoded)


def satellite_pass(Earthsatellite, start_date = datetime.now(), period_offset_days = 1, coords = wgs84.latlon(43.637 * N, 3.843 * E), elev_min = 0):
    ts = load.timescale()
    date = start_date
    t0 = ts.from_datetime(timezone('Europe/Paris').localize(date))  #- timedelta(hours=1)
    t1 = ts.from_datetime(timezone('Europe/Paris').localize(date + timedelta(days=period_offset_days)))
    t, events = Earthsatellite.find_events(coords, t0, t1, altitude_degrees=elev_min)
    event_names = [f'rise above {elev_min}', 'culminate', f'set below {elev_min}']
    eph = load('de421.bsp')
    sunlit = Earthsatellite.at(t).is_sunlit(eph)
    columns_head = ['Date(FR)', 'Event', 'Elevation', 'Illumination']
    Timestamp = []
    Elev = []
    Event = []
    Illu = []
    for ti, event, sunlit_flag in zip(t, events, sunlit):
        name = event_names[event]
        if sunlit_flag : state = 'in sunlight'
        else : state = 'in shadow'  
        difference = Earthsatellite - coords
        topocentric = difference.at(ti)
        alt, az, distance = topocentric.altaz()
        Timestamp.append((ti.astimezone(timezone('Europe/Paris'))).strftime('%d %b %Y - %H:%M:%S'))
        Elev.append(round(alt.degrees,1))
        Event.append(name)
        Illu.append(state)
    data_raw = [Timestamp, Event, Elev, Illu]
    dic_export = dict(zip(columns_head, data_raw))
    df = pd.DataFrame.from_dict(dic_export)
    df['Date(FR)'] = pd.to_datetime(df['Date(FR)'], dayfirst = True)
    return(df.set_index('Date(FR)'))

def est_visible(tle_line1 : str, tle_line2 : str, lat : float, lon : float, alt : float, date : datetime):
    ts = load.timescale()
    satellite = TLE_to_sat('\n'.join([tle_line1, tle_line2]))
    t = ts.utc(date.year, date.month, date.day, date.hour, date.minute, date.second)

    observateur = Topos(latitude_degrees=lat, longitude_degrees=lon, elevation_m=alt)

    # Calcul de la position
    difference = satellite - observateur
    topocentric = difference.at(t)
    alt, az, distance = topocentric.altaz()

    return alt.degrees

############################# DATE FORMAT #############################

def from_CEST_to_UTC(date):
    from pytz import timezone
    return(timezone('Europe/Paris').localize(date)).astimezone(tz=timezone('UTC'))

def from_UTC_to_CEST(date):
    from pytz import timezone
    return(timezone('UTC').localize(date)).astimezone(tz=timezone('Europe/Paris'))

############################## TLE GENERATION ##############################


def mean_motion_from_sma(sma_km):
    """Calculer la Mean Motion (révolution par jour) à partir du Semi-Major Axis (km)"""
    # Mean Motion en rad/s
    mean_motion_rad_s = math.sqrt(mu / (sma_km**3))

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
