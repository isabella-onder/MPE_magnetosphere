from astropy.io import fits
import numpy as np                                                                                                      
import matplotlib.pyplot as plt     


def many(eROdays, threshold, binsize):                                                                                                          
    file_list = ["/data36s/bella/erodat/erosita/eROday/" + "e_6_" + str(eROday) + "_002_c030.fits" for eROday in eROdays]
    fig, axs = plt.subplots(len(file_list)-1,2, sharex = 'col', sharey= True)
    for k in range(len(file_list) - 1):
        with fits.open(file_list[k]) as hdul_1, fits.open(file_list[k+1]) as hdul_2:

            bin_time = binsize

            data_1 = hdul_1[1].data 
            filtered_events_1 = data_1[data_1['PI'] > threshold]
            time_1 = filtered_events_1['TIME']
            reset_time_1 = time_1 - min(time_1)
            poubelles_1 = int((max(time_1)-min(time_1))/bin_time)


            data_2 = hdul_2[1].data 
            filtered_events_2 = data_2[data_2['PI'] > threshold]
            time_2 = filtered_events_2['TIME']
            reset_time_2 = time_2 - min(time_2)
            poubelles_2 = int((max(time_2)-min(time_2))/bin_time)

            #print('these are the poubelles_s', poubelles_1, poubelles_2)
            #print('these are the max(reset_time_s)', max(reset_time_1), max(reset_time_2))
            #print('these are the min(reset_time_s)', min(reset_time_1), min(reset_time_2))

            #print('these are the len(time_s)',len(time_1), len(time_2))

            n_1,bins_1 = np.histogram(reset_time_1, poubelles_1, (0,max(reset_time_1)))
            n_2,bins_2 = np.histogram(reset_time_2, poubelles_2, (0,max(reset_time_2)))
            #print('these are the len(n_s)',len(n_1), len(n_2))


            if len(n_2)!=len(n_1):
                print('eROdays' + str(eROdays[k])+ ' and ' + str(eROdays[k+1]) + ' do not have the same length')
                break

            n = n_2 - n_1
            normalised_n = n/bin_time

            #bins = np.delete(bins_1,0)
            time_in_hours = bins_1/3600

            bin_nb_for_histogram = 12

            axs[k][0].axhline(y=0,color ='grey', linestyle ='--')
            #axs[k][0].plot(time_in_hours, normalised_n, linewidth = 1)
            axs[k][0].stairs(normalised_n,time_in_hours)
            axs[k][1].axhline(y=0, color='grey', linestyle ='--')
            axs[k][1].hist(normalised_n, bins = bin_nb_for_histogram, orientation = 'horizontal')
    fig.suptitle('eRodays ' + str(eROdays[0]) + ' to ' + str(eROdays[-1]) + ' ; threshold (eV) = ' + str(threshold) + ' ; bins (s) = ' + str(binsize))
    fig.supxlabel("Time (hrs)")
    fig.supylabel("Counts/sec")
    plt.show()
