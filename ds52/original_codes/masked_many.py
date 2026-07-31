from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import search_around_sky
import astropy.units as u
import numpy as np
import scipy
import scipy.stats
import matplotlib.pyplot as plt





def many(eROdays, threshold, binsize):
    eROdays = list(range(eROdays[0], eROdays[-1]+1))
    file_list = ["/data36s/bella/erodat/erosita/eROday/" + "e_6_" + str(eROday) + "_002_c030.fits" for eROday in eROdays]
    fig, axs = plt.subplots(len(file_list)-1,2, sharex = 'col', sharey= True)
    for k in range(len(file_list) - 1):
        with fits.open(file_list[k]) as hdul_ero_1, fits.open(file_list[k+1]) as hdul_ero_2, fits.open("/data36s/bella/erodat/erosita/artxc_catalogue_ecliptic.fits") as hdul_cat:
            bin_time = binsize
            catalogue = hdul_cat[1].data
            mask_radius = 0.067*u.deg


            non_filtered_events_1 = hdul_ero_1[1].data
            events_1 = non_filtered_events_1[non_filtered_events_1['PI'] > threshold]
            events_elon_1 = events_1['ELON']
            events_elat_1 = events_1['ELAT']

            if events_elon_1[10000] < 180:
                min_threshold_1 = events_elon_1[10000] - 30
                max_threshold_1 = events_elon_1[10000] + 180 + 30

            elif events_elon_1[10000] > 180:
                min_threshold_1 = events_elon_1[10000] - 180 - 30
                max_threshold_1 = events_elon_1[10000] + 30

            #print(min_threshold_1, max_threshold_1) #the threshold makers seem to work pretty well
            cat_range_1 = catalogue[(catalogue['ELON'] > min_threshold_1) & (catalogue['ELON'] < max_threshold_1)]
            cat_elon_1 = cat_range_1['ELON']
            cat_elat_1 = cat_range_1['ELAT']

            events_coords_1 = SkyCoord(events_elon_1 * u.deg, events_elat_1 * u.deg, frame = 'barycentrictrueecliptic')
            cat_coords_1 = SkyCoord(cat_elon_1 * u.deg, cat_elat_1 * u.deg, frame = 'barycentrictrueecliptic')
            
            idx_cats_1, idx_counts_1,_ ,_ = search_around_sky(cat_coords_1,events_coords_1,mask_radius)
            Sco_X1_coords = SkyCoord(ra=[244.979]*u.degree, dec=[-15.640]*u.degree, frame ='icrs')
            mask_radius_Sco_X1 = 0.67*u.deg
            _, idx_counts_Sco_X1_1,_,_ = search_around_sky(Sco_X1_coords, events_coords_1,mask_radius_Sco_X1)

            idx_counts_1 = np.concatenate((idx_counts_1,idx_counts_Sco_X1_1))

            mask = np.ones(events_1.shape[0], dtype=bool)
            mask[np.unique(idx_counts_1)] = False
            filtered_events_1 = events_1[mask]

            time_1 = filtered_events_1['TIME']
            reset_time_1 = time_1 - min(time_1)
            poubelles_1 = int((max(time_1)-min(time_1))/bin_time)
            n_1,bins_1 = np.histogram(reset_time_1, poubelles_1, (min(reset_time_1),max(reset_time_1)))



            non_filtered_events_2 = hdul_ero_2[1].data
            events_2 = non_filtered_events_2[non_filtered_events_2['PI'] > threshold]
            events_elon_2 = events_2['ELON']
            events_elat_2 = events_2['ELAT']

            if events_elon_2[10000] < 180:
                min_threshold_2 = events_elon_2[10000] - 30
                max_threshold_2 = events_elon_2[10000] + 180 + 30

            elif events_elon_2[10000] > 180:
                min_threshold_2 = events_elon_2[10000] - 180 - 30
                max_threshold_2 = events_elon_2[10000] + 30

            #print(min_threshold_2, max_threshold_2) #the threshold makers seem to work pretty well
            cat_range_2 = catalogue[(catalogue['ELON'] > min_threshold_2) & (catalogue['ELON'] < max_threshold_2)]
            cat_elon_2 = cat_range_2['ELON']
            cat_elat_2 = cat_range_2['ELAT']

            events_coords_2 = SkyCoord(events_elon_2 * u.deg, events_elat_2 * u.deg, frame = 'barycentrictrueecliptic')
            cat_coords_2 = SkyCoord(cat_elon_2 * u.deg, cat_elat_2 * u.deg, frame = 'barycentrictrueecliptic')
            
            idx_cats_2, idx_counts_2,_ ,_ = search_around_sky(cat_coords_2,events_coords_2,mask_radius)
            Sco_X1_coords = SkyCoord(ra=[244.979]*u.degree, dec=[-15.640]*u.degree, frame ='icrs')
            mask_radius_Sco_X1 = 0.67*u.deg
            _, idx_counts_Sco_X1_2,_,_ = search_around_sky(Sco_X1_coords, events_coords_2,mask_radius_Sco_X1)

            idx_counts_2 = np.concatenate((idx_counts_2,idx_counts_Sco_X1_2))

            mask = np.ones(events_2.shape[0], dtype=bool)
            mask[np.unique(idx_counts_2)] = False
            filtered_events_2 = events_2[mask]

            time_2 = filtered_events_2['TIME']
            reset_time_2 = time_2 - min(time_2)
            poubelles_2 = int((max(time_2)-min(time_2))/bin_time)
            n_2,bins_2 = np.histogram(reset_time_2, poubelles_2, (min(reset_time_2),max(reset_time_2)))


            if len(n_2)!=len(n_1):
                print('eROdays' + str(eROdays[k])+ ' and ' + str(eROdays[k+1]) + ' do not have the same length')
                break

            '''Try to normalise before substracting instead, to do operations individually
            n = n_2 - n_1
            normalised_n = n/bin_time'''

            normalised_n_1 = n_1/bin_time
            normalised_n_2 = n_2/bin_time
            #print(normalised_n_1, normalised_n_2)
            normalised_n = normalised_n_2 - normalised_n_1
            
            idx_zeroes_1 = np.where(normalised_n_1 == 0)[0]
            #print('these are idx_zeroes_1',idx_zeroes_1)
            idx_zeroes_and_around_1 = list(idx_zeroes_1)
            for index in idx_zeroes_1:
                idx_zeroes_and_around_1.append(index-1)
                idx_zeroes_and_around_1.append(index+1)
            idx_zeroes_2 = np.where(normalised_n_2 == 0)[0]
            #print('these are idx_zeroes_2',idx_zeroes_2)
            idx_zeroes_and_around_2 = list(idx_zeroes_2)
            for index in idx_zeroes_2:
                idx_zeroes_and_around_2.append(index-1)
                idx_zeroes_and_around_2.append(index+1)
            idx_zeroes_and_around = idx_zeroes_and_around_1 + idx_zeroes_and_around_2
            idx_zeroes_and_around = list(set(idx_zeroes_and_around))
            #print('these are the lists combined', idx_zeroes_and_around)

            normalised_n_1_without_zeroes = np.delete(normalised_n_1,idx_zeroes_and_around)
            normalised_n_2_without_zeroes = np.delete(normalised_n_2,idx_zeroes_and_around)
            #normalised_n_without_zeroes_test = normalised_n_2_without_zeroes - normalised_n_1_without_zeroes 
            #print(normalised_n_1_without_zeroes, normalised_n_2_without_zeroes)
            
            
            time_in_hours = bins_1/3600

            bin_nb_for_histogram = 12
            
            '''#These are the statistical parameters WITH the zeroes
            mean = np.mean(normalised_n)
            median = np.median(normalised_n)
            mode = scipy.stats.mode(normalised_n)
            skew = scipy.stats.skew(normalised_n)
            print(mean, median, mode[0][0],skew)'''


            
            axs[k][0].axhline(y=0,color ='grey', linestyle ='--')
            axs[k][0].stairs(normalised_n,time_in_hours)
            if len(idx_zeroes_and_around) != 0:
                idx_zeroes_and_around_plus_one = idx_zeroes_and_around + [idx_zeroes_and_around[-1]+1]
                axs[k][0].plot(time_in_hours[idx_zeroes_and_around_plus_one], [0]*len(idx_zeroes_and_around_plus_one), 'r')


            normalised_n_without_zeroes = np.delete(normalised_n,idx_zeroes_and_around)
            #print(normalised_n_without_zeroes == normalised_n_without_zeroes_test)

            axs[k][1].axhline(y=0, color='grey', linestyle ='--')
            axs[k][1].hist(normalised_n_without_zeroes, bins = bin_nb_for_histogram, orientation = 'horizontal')

            #These are the statistical parameters WITHOUT zeroes
            mean = np.mean(normalised_n_without_zeroes)
            median = np.median(normalised_n_without_zeroes)
            mode = scipy.stats.mode(normalised_n_without_zeroes)
            skew = scipy.stats.skew(normalised_n_without_zeroes)
            print(mean, median, mode[0][0],skew)

            #axs[k][1].text(0.6,0.1, 'μ= ; med=  \n mode= ; skew= ', transform=axs[k][1].transAxes, size ='x-small')
    fig.suptitle('eRodays ' + str(eROdays[0]) + ' to ' + str(eROdays[-1]) + ' ; threshold (eV) = ' + str(threshold) + ' ; bins (s) = ' + str(binsize))
    fig.supxlabel('Time(h)          |          Counts')
    fig.supylabel('Counts/sec')
    #axs[4][0].set_title('Counts/sec')
    #axs[4][1].set_title('Counts')


    plt.show()
