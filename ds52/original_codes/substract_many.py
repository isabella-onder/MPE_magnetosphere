from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

def many(eROdays):
    file_list = ["/data36s/bella/erodat/erosita/eROday/" + "e_6_" + str(eROday) + "_002_c030.fits" for eROday in eROdays]
    fig, axs = plt.subplots(len(file_list)-1, sharex = True)
    for k in range(len(file_list) - 1):
        with fits.open(file_list[k]) as hdul_1, fits.open(file_list[k+1]) as hdul_2:
            bin_time = 32
            time_1 = hdul_1[1].data['TIME']
            reset_time_1 = time_1 - min(time_1)
            poubelles_1 = int((max(time_1)-min(time_1))/bin_time)
            

            time_2 = hdul_2[1].data['TIME']
            reset_time_2 = time_2 - min(time_2)
            poubelles_2 = int((max(time_2)-min(time_2))/bin_time)

            #print('these are the poubelles_s', poubelles_1, poubelles_2)
            #print('these are the max(reset_time_s)', max(reset_time_1), max(reset_time_2))
            #print('these are the min(reset_time_s)', min(reset_time_1), min(reset_time_2))

            #print('these are the len(time_s)',len(time_1), len(time_2))

            

            n_1,bins_1 = np.histogram(reset_time_1, poubelles_1, (0,max(reset_time_1)))
            n_2,bins_2 = np.histogram(reset_time_2, poubelles_2, (0,max(reset_time_2)))
            #print('these are the len(n_s)',len(n_1), len(n_2))
            

            n = n_2 - n_1
            normalised_n = n/bin_time

            bins = np.delete(bins_1,0)
            time_in_hours = bins/3600
            axs[k].plot(time_in_hours, normalised_n, linewidth = 1)
    fig.suptitle('eRodays ' + str(eROdays[0]) + ' to ' + str(eROdays[-1]))
    fig.supxlabel("Time (hrs)")
    fig.supylabel("Counts/sec")
    plt.show()
