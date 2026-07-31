import glob
import matplotlib.pyplot as plt
from astropy.io import fits


def obs_plotter(rev, obs_id):
    files = glob.glob(f"/xmm/archive/products/{rev}/{obs_id}/PN*/P*PN*FBKTSR*.FIT*")
    for filepath in files:
        with fits.open(filepath) as hdul:
            time = hdul[1].data['TIME']
            rate = hdul[1].data['RATE']

            print('the time bin is', time[10]-time[9], 'seconds long')

            fig, ax = plt.subplots()
            ax.set_xlabel("Time (XMM seconds)")
            ax.set_ylabel("Counts/sec XMM")
            ax.plot(time, rate, marker = '.')

            plt.show()

obs_plotter('3641', '0823810201')
