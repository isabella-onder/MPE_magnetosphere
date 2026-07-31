from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import search_around_sky
from astropy.stats import bayesian_blocks
import astropy.stats
import astropy.units as u
import numpy as np
import scipy
import scipy.stats
import matplotlib.pyplot as plt
import json

with fits.open('/data36s/bella/erodat/xmm/xsa_archive_2019_2022_repaired.fits') as hdul:
    data = hdul[1].data

    obsIDs = data['OBSERVATION']
    starts = data['TIME_0_XMM']
    ends = data['TIME_E_XMM']
    
    output = []
    for start, end, obs in zip(starts, ends, obsIDs):
        try:
            duration = int(float(end)) - int(float(start)) 
        except ValueError:
            continue
        overarching = [obs, duration] 
        output.append(overarching)

    with open("obsID_durations.json", "w") as f:
        json.dump(output, f, indent = 2)
        print('the list containing [[obsID, duration], ...] has been dumped into obsID_durations.json')

