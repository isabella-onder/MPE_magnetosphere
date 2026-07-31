from astropy.io import fits
import numpy as np
from datetime import datetime, timedelta

global starting_point 
starting_point = datetime(1999,12,31,21,00,00)

def time_ETC_to_eRODAY(time_eROSITA_ETC):
    time_eROSITA_ETC = datetime.fromisoformat(time_eROSITA_ETC)
    print(time_eROSITA_ETC)
    hours_passed = (time_eROSITA_ETC - starting_point).total_seconds() / 3600 + 4
    eroday = hours_passed / 4
    return eroday

with fits.open("/data36s/bella/erodat/erosita/P_PLAN_628608605_631297805_HB.fits") as hdul:
    UTC_times_start =hdul[1].data['TSTART_UTC']
    UTC_times_end = hdul[1].data['TSTOP_UTC']
    #print(UTC_times)
    eROdays = []
    for k in range(len(UTC_times_start)):
        eROday_start = int(time_ETC_to_eRODAY(UTC_times_start[k]))
        eROday_end = int(time_ETC_to_eRODAY(UTC_times_end[k]))
        eROday_range = list(range(eROday_start, eROday_end + 1))
        for eROday in eROday_range:
            eROdays.append(eROday)
    eROdays = list(set(eROdays))
    print(eROdays, len(eROdays))


with fits.open("/data36s/bella/erodat/erosita/erosita_artxc_fieldscans.fits") as hdul:
    eROdays_beginning = hdul[1].data['eROday_b']
    eROdays_end = hdul[1].data['eROday_e']
    #extracting the erodays from the artxc_field_scans file, to then include in the invalid list for detector
    eROdays_fieldscan = []
    for k in range (len(eROdays_end)):
        fieldscan = list(range(int(eROdays_beginning[k]), int(eROdays_end[k]) + 1))
        eROdays_fieldscan.extend(fieldscan)
    fieldscan_days = list(set(eROdays_fieldscan))
    print('fieldscan days',fieldscan_days)

