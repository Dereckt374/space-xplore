import os

from pytz import timezone
from skyfield.api import EarthSatellite, load, wgs84, N, E
from datetime import datetime, timedelta, time
import requests
import pandas as pd
import sys
import requests
import pandas as pd

path = '/home/virgil/Documents/Git/gseg/tools/functions'
sys.path.append(path)

from fct_orbit import *

MCC_URL = "http://162.38.203.15:8000"

date = datetime.now()
day = date.day
minutes = date.minute
ts = load.timescale()
t = ts.from_datetime(timezone("Europe/Paris").localize(date))
t0 = datetime.now()
offset_day = 2
minimal_elev = 5
loc = wgs84.latlon(43.637 * N, 3.843 * E)

def transf(df): 
    lst = []
    for index, row in df.iterrows():
        name = str(row['Satellite']) + ' ' + str(row['Elev.']) + '° - ' + index
        end = (datetime.strptime(index, "%H:%M") + timedelta(minutes = 20)).strftime("%H:%M")
        date = f'- {index} - {end} '
        lst.append(date + name)
    return(lst)

if __name__ == "__main__":
    all_passes = get_operable_passes(MCC_URL, t0, loc, offset_day, minimal_elev)
    for i in transf(all_passes):
        print(i)

