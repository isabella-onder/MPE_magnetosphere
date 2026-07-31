from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import search_around_sky
import astropy.units as u
import numpy as np
import scipy
import scipy.stats
import matplotlib.pyplot as plt

def lightcurve_plotter(eROday, threshold, binsize):
    filename = "e_6_" + str(eROday) + "_002_c030.fits"
    filepath = "/data36s/bella/erodat/erosita/eROday/" + filename
    fig,axs = plt.subplots(1,2,sharey=True)

    #threshold = 4500

    with fits.open(filepath) as hdul_ero, fits.open("/data36s/bella/erodat/erosita/artxc_catalogue_ecliptic.fits") as hdul_cat:
        #events = hdul_ero[1].data
        non_filtered_events = hdul_ero[1].data
        events = non_filtered_events[non_filtered_events['PI'] > threshold]
        catalogue = hdul_cat[1].data

        events_elon = events['ELON']
        events_elat = events['ELAT']

        '''if events_elon[10000] < 180:
            min_threshold = events_elon[10000] - 30
            max_threshold = events_elon[10000] + 180 + 30

        elif events_elon[10000] > 180:
            min_threshold = events_elon[10000] - 180 - 30
            max_threshold = events_elon[10000] + 30

        print(min_threshold, max_threshold) #the threshold makers seem to work pretty well
        cat_range = catalogue[(catalogue['ELON'] > min_threshold) & (catalogue['ELON'] < max_threshold)]
        cat_elon = cat_range['ELON']
        cat_elat = cat_range['ELAT']'''

        mean = np.mean(events_elon)
        min_threshold = mean - 120
        max_threshold = mean + 120


        cat_range = catalogue[((catalogue['ELON'] > min_threshold) & (catalogue['ELON'] < min_threshold + 60)) |((catalogue['ELON'] > max_threshold - 60) & (catalogue['ELON'] < max_threshold))]
        cat_elon = cat_range['ELON']
        cat_elat = cat_range['ELAT']

        events_coords = SkyCoord(events_elon * u.deg, events_elat * u.deg, frame = 'barycentrictrueecliptic')
        cat_coords = SkyCoord(cat_elon * u.deg, cat_elat * u.deg, frame = 'barycentrictrueecliptic')
        #print('these are the cat_coords (elon,elat)', cat_elon, cat_elat)
        #print('these are the event_coords (elon, elat)',events_elon, events_elat)
        mask_radius = 0.067*u.deg
        idx_cats, idx_counts,_ ,_ = search_around_sky(cat_coords,events_coords,mask_radius)
        #print('these are the mapping indexes lists from search_around_sky',idx_cats,idx_counts, 'and their lengths', len(idx_cats), len(idx_counts))
        
        #for Sco-X1 specifically
        Sco_X1_coords = SkyCoord(ra=[244.979]*u.degree,dec=[-15.640]*u.degree, frame ='icrs')
        mask_radius_Sco_X1 = 3*u.deg
        _, idx_counts_Sco_X1,_,_ = search_around_sky(Sco_X1_coords, events_coords,mask_radius_Sco_X1)
        
        idx_counts = np.concatenate((idx_counts,idx_counts_Sco_X1))





        mask = np.ones(events.shape[0], dtype=bool)
        mask[np.unique(idx_counts)] = False
        filtered_events = events[mask]
        
        print(len(events), len(filtered_events))

        time = filtered_events['TIME']
        reset_time = time - min(time)
        bin_time = binsize
        poubelles = int((max(time)-min(time))/bin_time)
        print('these are poubelles', poubelles)
        n,bins = np.histogram(reset_time, poubelles, (min(reset_time),max(reset_time)))
        normalised_n = n/bin_time
        print(normalised_n)
        print('these are my lengths for n and bins', len(n), len(bins), 'these are my normalised_n', normalised_n, 'these are my bins')
        #print(normalised_n)

        #idx_zeroes = np.where(normalised_n == 0)[0]
        #normalised_n = np.delete(normalised_n,idx_zeroes)
        
        time_in_hours = bins/3600

        bin_nb_for_histogram = 12

        ''' Parameters WITH the zeroes
        mean = np.mean(normalised_n)
        median = np.median(normalised_n)
        mode = scipy.stats.mode(normalised_n)
        skew = scipy.stats.skew(normalised_n)
        print('Mean: ',mean,'\nMedian: ',  median,'\nMode: ', mode[0][0],'\nSkew: ',skew)'''

        axs[0].stairs(normalised_n, time_in_hours)

        #NB: even though the zeroes still appear on the stairs plot, they are NOT included in the histogram and statistical parameters' calculationsi
        idx_zeroes = np.where(normalised_n == 0)[0]
        idx_zeroes_and_around = list(idx_zeroes)
        for index in idx_zeroes:
            idx_zeroes_and_around.append(index-1)
            idx_zeroes_and_around.append(index+1)
        idx_zeroes_and_around = list(set(idx_zeroes_and_around))


        normalised_n_without_zeroes = np.delete(normalised_n,idx_zeroes_and_around)

        #parameters WITHOUT the zeroes
        mean= np.mean(normalised_n_without_zeroes)
        median = np.median(normalised_n_without_zeroes)
        mode = scipy.stats.mode(normalised_n_without_zeroes)
        skew = scipy.stats.skew(normalised_n_without_zeroes)
        print('Mean: ',mean,'\nMedian: ',  median,'\nMode: ', mode[0][0],'\nSkew: ',skew)

        if len(idx_zeroes_and_around)!=0:
            idx_zeroes_and_around_plus_one = idx_zeroes_and_around + [idx_zeroes_and_around[-1]+1]
            axs[0].plot(time_in_hours[idx_zeroes_and_around_plus_one], [mean]*len(idx_zeroes_and_around_plus_one), 'r')
        axs[1].hist(normalised_n_without_zeroes, bins = bin_nb_for_histogram, orientation ='horizontal')
        #axs[1].text(0.6,0.1, 'μ='+str(mean)+' ; med='+str(median)+'  \n mode='+str(mode)+ ' ; skew='+str(skew), transform=axs[1].transAxes, size ='x-small')
        test_1, test_2 = np.histogram(normalised_n_without_zeroes, bins = bin_nb_for_histogram)
        
       # print('test 1', test_1, ' \n test 2', test_2)

        #plt.plot(time_in_hours, normalised_n, linewidth = 1)
        fig.supxlabel("Time(h)")
        fig.supylabel("Counts/sec")
        fig.suptitle('eROday ' +str(eROday)+ ' ; threshold (eV) =  ' +str(threshold)+ ' ; binsize (s) = ' +str(binsize)+ ' with mask radius (degrees) = ' +str(mask_radius))
        plt.show()
