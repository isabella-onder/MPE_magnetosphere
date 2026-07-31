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
    with fits.open(filepath) as hdul_ero, fits.open("/data36s/bella/erodat/erosita/artxc_catalogue_ecliptic.fits") as hdul_art, fits.open("/data36s/bella/erodat/erosita/maxi_bright_src.fits") as hdul_maxi:
        non_filtered_events = hdul_ero[1].data                                                                          
        events = non_filtered_events[non_filtered_events['PI'] > threshold]   
        

        if len(events)==0:
            print('Len(events)==0. Something has gone very wrong.')
            mean, median, mode, skew = -1000, -1000, -1000, -1000 #putting 1000 as an absurd value to flag it and leave it
            return mean, median, mode, skew


        catalogue_art = hdul_art[1].data
        catalogue_maxi = hdul_maxi[1].data

        events_elon = events['ELON']                                                                                           
        events_elat = events['ELAT']

        mean = np.mean(events_elon)
        min_threshold = mean - 120
        max_threshold = mean + 120

        cat_range_art = catalogue_art[((catalogue_art['ELON'] > min_threshold) & (catalogue_art['ELON'] < min_threshold + 60)) |((catalogue_art['ELON'] > max_threshold - 60) & (catalogue_art['ELON'] < max_threshold))]
        cat_elon_art = cat_range_art['ELON']
        cat_elat_art = cat_range_art['ELAT']

        cat_range_maxi = catalogue_maxi[((catalogue_maxi['LAMDA'] > min_threshold) & (catalogue_maxi['LAMDA'] < min_threshold + 60)) |((catalogue_maxi['LAMDA'] > max_threshold - 60) & (catalogue_maxi['LAMDA'] < max_threshold))]
        cat_elon_maxi = cat_range_maxi['LAMDA']
        cat_elat_maxi = cat_range_maxi['BETA']
        

        #STRONG SUSPICION HE SWAPPED LAMDA AND BETA
        #cat_coords_test = SkyCoord(cat_elon_maxi * u.deg, cat_elat_maxi * u.deg, frame = 'barycentrictrueecliptic')

        cat_elon = np.concatenate((cat_elon_art,cat_elon_maxi))
        cat_elat = np.concatenate((cat_elat_art, cat_elat_maxi))

        events_coords = SkyCoord(events_elon * u.deg, events_elat * u.deg, frame = 'barycentrictrueecliptic')
        cat_coords = SkyCoord(cat_elon * u.deg, cat_elat * u.deg, frame = 'barycentrictrueecliptic')

        mask_radius = 0.067*u.deg
        idx_cats, idx_counts,_ ,_ = search_around_sky(cat_coords,events_coords,mask_radius)

        #for Sco-X1 specifically
        Sco_X1_coords = SkyCoord(ra=[244.979]*u.degree,dec=[-15.640]*u.degree, frame ='icrs')
        mask_radius_Sco_X1 = 3*u.deg
        _, idx_counts_Sco_X1,_,_ = search_around_sky(Sco_X1_coords, events_coords,mask_radius_Sco_X1)
        idx_counts = np.concatenate((idx_counts,idx_counts_Sco_X1))

        mask = np.ones(events.shape[0], dtype=bool)
        mask[np.unique(idx_counts)] = False
        filtered_events = events[mask]

        time = filtered_events['TIME']
        reset_time = time - min(time)
        bin_time = binsize
        poubelles = int((max(time)-min(time))/bin_time)
        n,bins = np.histogram(reset_time, poubelles, (min(reset_time),max(reset_time)))
        normalised_n = n/bin_time

        time_in_hours = bins/3600

        bin_nb_for_histogram = 12

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
        return(mean,median,mode,skew, normalised_n_without_zeroes)
        print('Mean: ',mean,'\nMedian: ',  median,'\nMode: ', mode[0][0],'\nSkew: ',skew)                                        
            


def many(range_eROday, threshold, binsize): #range is a list with starting and finishing eROday, threshold in eV and binsize in secs
    eROday_range_list = list(range(range_eROday[0], range_eROday[-1]+1))
    file_list = ["/data36s/bella/erodat/erosita/eROday/" + "e_6_" + str(eROday) + "_002_c030.fits" for eROday in eROday_range_list]    
    means = []
    medians = []
    modes = []
    skews = []                               
    for k in range(len(file_list) - 1):                                             
        with fits.open(file_list[k]) as hdul_ero_1, fits.open(file_list[k+1]) as hdul_ero_2, fits.open("/data36s/bella/erodat/erosita/artxc_catalogue_ecliptic.fits") as hdul_cat:  
            bin_time = binsize                                                          
            catalogue = hdul_cat[1].data                                                
            mask_radius = 0.067*u.deg

            non_filtered_events_1 = hdul_ero_1[1].data
            events_1 = non_filtered_events_1[non_filtered_events_1['PI'] > threshold]
            events_elon_1 = events_1['ELON']
            events_elat_1 = events_1['ELAT']

            mean_1 = np.mean(events_elon_1)
            min_threshold_1 = mean_1 - 120
            max_threshold_1 = mean_1 + 120

            cat_range_1 = catalogue[((catalogue['ELON'] > min_threshold_1) & (catalogue['ELON'] < min_threshold_1 + 60)) |((catalogue['ELON'] > max_threshold_1 - 60) & (catalogue['ELON'] < max_threshold_1))]
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

            mean_2 = np.mean(events_elon_2)
            min_threshold_2 = mean_2 - 120
            max_threshold_2 = mean_2 + 120

            cat_range_2 = catalogue[((catalogue['ELON'] > min_threshold_2) & (catalogue['ELON'] < min_threshold_2 + 60)) |((catalogue['ELON'] > max_threshold_2 - 60) & (catalogue['ELON'] < max_threshold_2))]
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
            normalised_n_1_without_zeroes = np.delete(normalised_n_1,idx_zeroes_and_around)
            normalised_n_2_without_zeroes = np.delete(normalised_n_2,idx_zeroes_and_around)

            time_in_hours = bins_1/3600

            bin_nb_for_histogram = 12


    

            normalised_n_without_zeroes = np.delete(normalised_n,idx_zeroes_and_around)
        
            #These are the statistical parameters WITHOUT zeroes
            mean = np.mean(normalised_n_without_zeroes)
            median = np.median(normalised_n_without_zeroes)
            mode = scipy.stats.mode(normalised_n_without_zeroes)
            skew = scipy.stats.skew(normalised_n_without_zeroes)
        
            means.append(mean)
            medians.append(median)
            modes.append(mode)
            skews.append(skew)
    return(means,medians,modes,skews)

            #might have to make a list for it to return in

    

def detect(range_of_eROdays, threshold, bins):
    eROday_list = list(range(range_of_eROdays[0], range_of_eROdays[-1]+1))

    single_means = []
    single_medians = []
    single_modes = []
    single_skews = []

    mean_threshold = 2.0
    skew_threshold = 0.9

    enhancements = []  #for the time being, we search for large, notable enhancements
    swab = [] #need for further test: high mean, low skew, may just be background if Gaussian nonetheless

    for eROday in eROday_list:
        try:
            single_mean, single_median, single_mode, single_skew, normalised_n = lightcurve_plotter(eROday, threshold, bins)

        except FileNotFoundError:
            continue
        except ValueError:
            continue

        if single_mean == -1000:
            continue


        single_means.append(single_mean)
        single_medians.append(single_median)
        single_modes.append(single_mode[0][0])
        single_skews.append(single_skew)
        print(eROday, 'has been completed')

        if single_mean > mean_threshold:
            enhancements.append(eROday)
            print('eROday'+str(eROday)+' had mean cts/s > '+str(mean_threshold))
        
        if single_skew > skew_threshold:
            enhancements.append(eROday)
            print('eROday'+str(eROday)+' had skew > '+str(skew_threshold))
    #print(single_means, single_medians, single_modes, single_skews)


        if (single_skew < skew_threshold and single_mean > mean_threshold):
            swab.append(eROday)
    
    enhancements = list(set(enhancements))
    enhancements.sort()
    
    if len(enhancements)==0:
        print('No enhancements detected in this range')
    else:
        print('Potential list of enhancements: ',enhancements, '(nb = ',len(enhancements),')')
        print('Swab list (require further testing):', swab)
    return(enhancements, swab)



def mean_colorbar(range_of_eROdays,threshold,bins): #need to decide whether this is to see periods of enhancements or background
                                                    #for the former: would keep threshold high, to see visually enhancement periods
                                                    #for the later: keep threshold around cutoff for enhancements (e.g.2.1/2.2) to see how background changes

    eROday_list = list(range(range_of_eROdays[0], range_of_eROdays[-1]+1))

    means = []
    missing_eROdays = []
    issue = []
    suspicious_high = []
    suspicious_low = []
    bg_means = []
    potential_large_enhancements = []

    
    for eROday in eROday_list:
        try:
            single_mean, _, _, _, _ = lightcurve_plotter(eROday, threshold, bins)
            means.append(single_mean)
            print('eROday',eROday,'has been computed')
            if 2.2 <single_mean < 10:
                print('eROday',eROday,'has ENHANCEMENT LIKE mean',single_mean)
                potential_large_enhancements.append(eROday)
            elif  10 < single_mean < 999 :
                print('eROday', eROday, 'has SUSPICIOUSLY HIGH mean', single_mean)
                suspicious_high.append(eROday)
            elif single_mean == -1000 :
                print('eROday', eROday, 'has no events')
                issue.append(eROday)
            elif single_mean < 1.0:
                print('eROday', eROday, 'has SUSPICTIOUSLY LOW mean', single_mean)
                suspicious_low.append(eROday)
            else:
                bg_means.append(single_mean)
    


        except FileNotFoundError:
            means.append(0)
            missing_eROdays.append(eROday)
            continue
        except ValueError:
            means.append(0)
            issue.append(eROday)
            continue
    
    n = 30 #number of days in a month (i.e. by how much we want to group them up)
    print(len(means))
    grouped_means = [means[i:i + n] for i in range(0, len(means), n)]
    while len(grouped_means[-1]) != n:
        grouped_means[-1].append(0)

    #print(grouped_means)
    print('The following eROday files were missing', missing_eROdays)
    print('The following eROdays have issues', issue)
    print('The following eROdays are suspiciously high', suspicious_high)
    print('The following eROdays are suspiciously low', suspicious_low)
    print('The following eROdays potentially have large enhancements', potential_large_enhancements)

    fig, ax = plt.subplots()
    cmap = plt.cm.viridis
    cmap.set_under('black')
    cmap.set_over('black')
    im_hopefully = ax.imshow(grouped_means, vmin = min(bg_means), vmax = max(bg_means))
    plt.colorbar(im_hopefully)
    plt.title('Mean counts/s for every eROday')
    plt.ylabel('Rows = approx one month')
    plt.xlabel('Days')
    plt.show()

