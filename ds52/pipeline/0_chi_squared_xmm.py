from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import search_around_sky
import astropy.units as u
import numpy as np
import scipy
import scipy.stats
import matplotlib.pyplot as plt
import json

def lightcurve_plotter(eROday, threshold, binsize):
    filename = "e_6_" + str(eROday) + "_002_c030.fits"
    filepath = "/data36s/bella/erodat/erosita/eROday/" + filename                                                                                                                                                                                                                                                                                                                                                                
    with fits.open(filepath) as hdul_ero, fits.open("/data36s/bella/erodat/erosita/artxc_catalogue_ecliptic.fits") as hdul_art, fits.open("/data36s/bella/erodat/erosita/maxi_bright_src.fits") as hdul_maxi:
        non_filtered_events = hdul_ero[1].data                                                                          
        events = non_filtered_events[non_filtered_events['PI'] > threshold]   
    

        if len(events)==0:
            print('Len(events)==0. Something has gone very wrong with eROday'+str(eROday))
            mean, median, mode, skew = -1000, -1000, -1000, -1000 #putting 1000 as an absurd value to flag it and leave it
            return mean, median, mode, skew, -1000


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
        
        if len(normalised_n_without_zeroes)==0:
            print('Len(normalised_events_without_zeroes)==0. There are NO counts once stars are substracted - that is suspicious, check len(events) - eROday'+str(eROday))
            mean, median, mode, skew = -1000, -1000, -1000, -1000 #putting 1000 as an absurd value to flag it and leave it
            return mean, median, mode, skew, -1000

        max_reasonable_count = 15
        cutoff_test = all(x < max_reasonable_count for x in normalised_n_without_zeroes)
        if not cutoff_test:
            mean,median,mode,skew = -1000,-1000,-1000,-1000
            print('there is a normalised n without zeroes >', max_reasonable_count, ' eROday ', eROday, '  was suspicious and therefore removed')
            return mean,median,mode,skew,-1000

        
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
            
            print(len(n_1), len(n_2))
            if len(n_2)!=len(n_1):
                print('eROdays' + str(eROday_range_list[k])+ ' and ' + str(eROday_range_list[k+1]) + ' do not have the same length: ' +str(len(n_1)) + ' vs ' +str(len(n_2)) )
                means, medians, modes, skews = [1000],[1000],[1000],[1000] #putting 1000 as absurd values, to ensure detection 
                return (means, medians, modes, skews)

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

import scipy.stats as stats

def chi2(eROday, threshold, bins):
    plot = False

    _,_,_,_,data = lightcurve_plotter(eROday,threshold,bins)
    mu, std = np.mean(data), np.std(data)

    # Histogram (observed frequencies)
    observed, bin_edges = np.histogram(data,16)
    N = len(data)

    # Expected frequencies from fitted Gaussian
    expected_probs = stats.norm.cdf(bin_edges[1:], mu, std) - stats.norm.cdf(bin_edges[:-1], mu, std)
    expected = expected_probs * N

    # Chi-square statistic
    chi2 = np.sum((observed - expected) ** 2 / expected)
    dof = len(observed) - 2  # 2 fitted params: mu, std
    reduced_chi2 = chi2/dof
    p_value = stats.chi2.sf(chi2, dof)

    if plot:
        plt.figure(figsize=(6,4))
        plt.hist(data, bins=bin_edges, alpha=0.6, label="Observed", density=False)
        x = np.linspace(min(data), max(data), 200)
        pdf = stats.norm.pdf(x, mu, std) * N * (bin_edges[1]-bin_edges[0])
        plt.plot(x, pdf, 'r-', label=f"Gaussian fit\nμ={mu:.2f}, σ={std:.2f}")
        plt.xlabel("Counts/sec")
        plt.ylabel("Frequency")
        plt.legend()
        plt.show()
    print(eROday, chi2, reduced_chi2, dof, p_value)
    return chi2, reduced_chi2, dof, p_value

def consecutive(data, stepsize=1):
    return np.split(data, np.where(np.diff(data)!=stepsize)[0]+1)

def detect(range_of_eROdays, threshold, bins):
    eROday_list = list(range(range_of_eROdays[0], range_of_eROdays[-1]+1))
    

    filter_wheel_eROdays = list(range(44115,45081,42))
    fieldscan_eROdays = [43008, 43009, 43010, 43011, 43013, 43014, 43015, 43016, 43017, 43019, 43020, 43021, 43022, 43023, 43025, 43026, 43027, 43028, 43029, 43031, 43032, 43033, 43034, 43035, 43043, 43044, 43045, 43046, 43047, 43049, 43050, 43051, 43052, 43053, 43054, 43055, 43056, 43057, 43058, 43059, 43060, 43061, 43062, 43063, 43064, 43065, 43066, 43067, 43068, 43069, 43070, 43071, 43072, 43073, 43074, 43075, 43076, 43077, 43078, 43079, 43080, 43081, 43082, 43083, 43086, 43087, 43088, 43089, 43090, 43092, 43093, 43094, 43095, 43096, 43098, 43099, 43100, 43101, 43102, 43104, 43105, 43106, 43107, 43108, 43111, 43112, 43113, 43114, 43117, 43118, 43119, 43120, 43123, 43124, 43125, 43126, 43129, 43130, 43131, 43132, 43135, 43136, 43137, 43138, 43140, 43141, 43142, 43143, 43144, 43146, 43147, 43148, 43149, 43150, 43152, 43153, 43154, 43155, 43156, 43158, 43159, 43160, 43161, 43162, 43165, 43166, 43167, 43168, 43171, 43172, 43173, 43174, 43177, 43178, 43179, 43180, 43182, 43183, 43184, 43185, 43195, 43196, 43197, 43198, 43201, 43202, 43203, 43204, 43206, 43207, 43208, 43209, 43210, 43212, 43213, 43214, 43215, 43218, 43219, 43220, 43221, 43238, 43239, 43240, 43242, 43243, 43244, 43245, 43246, 43248, 43249, 43250, 43251, 43252, 43254, 43255, 43256, 43257, 43258, 43260, 43261, 43262, 43263, 43264, 43278, 43279, 43280, 43281, 43282, 43284, 43285, 43286, 43287, 43288, 43324, 43325, 43326, 43327, 43328, 43329, 43330, 43331, 43332, 43333, 43352, 43353, 43354, 43355, 43356, 43357, 43358, 43359, 43366, 43367, 43368, 43369, 43370, 43371, 43372, 43373, 43374, 43380, 43381, 43382, 43383, 43384, 43477, 43478, 43479, 43480, 43481, 43482, 43483, 43484, 43485, 43486, 43487, 43488, 43489, 43490, 43491, 43492, 43493, 43494, 43495, 43496, 43497, 43498, 43499, 43500, 43501, 43502, 43530, 43531, 43532, 43533, 43535, 43536, 43537, 43538, 43539, 43540, 43541, 43542, 43543, 43544, 43545, 43555, 43556, 43557, 43558, 43559, 43560, 43561, 43562, 43563, 43564, 43565, 43566, 43567, 43568, 43569, 43570, 43571, 43637, 43638, 43639, 43640, 43641, 43653, 43654, 43655, 43663, 43664, 43665, 43666, 43667, 43668, 43669, 43670, 43671, 43672, 43673, 43674, 43675, 43676, 43677, 43678, 43679, 43680, 43681, 43682, 43683, 45798, 45799, 45800, 45801, 45802, 44837, 44838, 44839, 44840, 44841, 44842, 44843, 44844, 44845, 44846, 44847, 44848, 44849, 42969, 42993]
    invalid_eROdays_set = set(filter_wheel_eROdays + fieldscan_eROdays) #to disregard all invalid eROdays (which we chose as all field scans etc.)
# wrote it this way so that it would be easy to add different lists of different invalid reasons

    single_means = []
    single_medians = []
    single_modes = []
    single_skews = []

    CME_mean_threshold = 2.5
    #mean_threshold = 2.0

    if bins == 32:
        skew_threshold = 0.4
        print('the bins are 32, the chosen threshold was 0.4')
    elif bins == 120:
        skew_threshold = 0.7
        print('the bins are 120, the chosen threshold was 0.7')
    else: #default
        skew_threshold = 0.7
        print('the default threshold was chosen')

    enhancements = []  #for the time being, we search for large, notable enhancements
    swab_high_mean = [] #need for further test: high mean, low skew, may just be background if Gaussian nonetheless
    skewed_but_average = [] #need for further test: peaks without high background. If consecutive and when substracted no skew, probably just source - otherwise genuine peak
    CME_or_lag = [] #to catch days with very high average counts: either very very active OR something particularly odd
    for eROday in eROday_list:
       # print('this', eROday,' is being processed')
        if eROday in invalid_eROdays_set:
           # print(eROday, 'was in invalid day')
            continue

                
        if (eROday > 43411 and eROday < 46309):         #accomodating for the different setups and therefore different background rates
           # print('eroday check', eROday, 'between 43411 and 46309')
            mean_threshold = 1.801
        else:
           # print('eroday check',eROday,'below or equal 43411 and above or equal 46309')
            mean_threshold = 2.001
        
        try:
           # print('try was tried')
            single_mean, single_median, single_mode, single_skew, normalised_n = lightcurve_plotter(eROday, threshold, bins)
        

        except FileNotFoundError:
            continue
        except ValueError:
            continue

        if single_mean == -1000:
            print(eROday, 'has single mean -1000')
            continue

        print(eROday, single_mean)
        single_means.append(single_mean)
        single_medians.append(single_median)
        single_modes.append(single_mode[0][0])
        single_skews.append(single_skew)

        if (single_mean > mean_threshold and single_skew > skew_threshold):
            enhancements.append(eROday)
            print('eROday'+str(eROday)+' had mean cts/s > '+str(mean_threshold)+' & skew > ' +str(skew_threshold))    

        if (single_skew < skew_threshold and single_mean > mean_threshold):
            swab_high_mean.append(eROday)
            print('eROday'+str(eROday)+' had mean cts/s > '+str(mean_threshold)+' & skew < ' +str(skew_threshold))
            
            try:
                chi_squared, chi_squared_nu, degrees_of_freedom, p_value = chi2(eROday, threshold, 32)

            except TypeError: #this is in case doing it with threshold 15 counts in smaller bins is passed: if any -1000 normalised_values is returned, then pass anyways (since there is some type of error and we want to disregard the eROday)
                continue

            if chi_squared_nu > 2.9:#because this would indicate it is NOT Gaussian. May need to change threshold/
                #print(eROday, chi_squared_nu, 'AAAAAAAAAAAAAAAAA')
                enhancements.append(eROday)

        if single_mean > CME_mean_threshold:
            CME_or_lag.append(eROday)

        if (single_skew > skew_threshold and single_mean < mean_threshold):
            skewed_but_average.append(eROday)
            print ('eROday'+str(eROday)+' had mean cts/s < '+str(mean_threshold)+' & skew > ' +str(skew_threshold))

        print(eROday,' has been completed')
    
    enhancements_av_skewed = []
    #skewed_but_average = [48509,48510,48511,48512,48513,48514,48515]
    if not skewed_but_average: #if there are no skewed_but_average days at all, 
        consecutives = []
    else:
        consecutives = consecutive(skewed_but_average)
    print(consecutives)
    
    #print(consecutives, len(consecutives))

    #split_up = consecutive(skewed_but_average)
    #consecutives = [number for number in split_up if len(number)>1]
    #singles = [number for number in split_up if len(number)==1]
    #enhancements_av_skewed.append(singles)
    
    if consecutives: #if there are no skewed but average days (here named consecutives, but they can be of length 1), code does not even try running the many code
        for k in range (len(consecutives)):  #to ensure that '1day' consecutives i.e. stand alones are kept
            if len(consecutives[k]) == 1:
                print('there was a len(1) consecutives indeed and its type is', type(consecutives[k]))
                single_consecutive = consecutives[k][0]
                enhancements_av_skewed.append(single_consecutive)
                continue
            means_c, medians_c, modes_c, skews_c = many([consecutives[k][0], consecutives[k][-1]], threshold, bins)
            if means_c == [1000]:
                continue
            skew_c_threshold = 0.5 
            print(skews_c,'skews_c')
            for i in range(len(skews_c)):
                if abs(skews_c[i]) > skew_c_threshold:
                    print('this is skews_c[i]', skews_c[i], 'at index i', i,'and consecutives[k][i+1] will be added', consecutives[k][i+1]) 
                    enhancements_av_skewed.append(consecutives[k][i+1]) # because in sub n_2 - n_1 we want n_2 if it is high skew

            #test_c = all(x < skew_c_threshold for x in skews_c)
            #if test_c == False:
                #enhancements_av_skewed.append(consecutives[k])
            #enhancements.append(eROday)
            #print('eROday'+str(eROday)+' had skew > '+str(skew_threshold))
    #print(single_means, single_medians, single_modes, single_skews)
    #print(enhancements_av_skewed, type(enhancements_av_skewed), 'AAAA')
    #enhancements_av_skewed_flat = [x for xs in enhancements_av_skewed for x in xs] #to flatten enhancements from the average and skewed
    enhancements = enhancements + enhancements_av_skewed
    #print(enhancements)
    enhancements = list(set(enhancements))
    enhancements.sort()
    enhancements = [int(x) for x in enhancements] #they are just eROdays so to make them json supscritable just take the integer
    
    
    if len(enhancements)==0:
        print('No enhancements detected in this range')
    else:
        print('Potential list of enhancements: ',enhancements, '(nb = ',len(enhancements),')', 'type: ', type(enhancements))
        print('Potential list of STRONG enhancements:', CME_or_lag, '(nb = ',len(CME_or_lag),')')
        print('Swab list (high mean & low skew, submitted to Gaussian test):', swab_high_mean)
        print('Sources or skewed flares (average mean & high skew, submitted to consecutive test):', skewed_but_average)
    with open("all_enhancements_output.json", "w") as f:
        json.dump(enhancements, f)
    with open("potential_strong_enhancements_output.json", "w") as f:
        json.dump(CME_or_lag,f)
    with open("high_mean_low_skew_output.json", "w") as f:
        json.dump(swab_high_mean, f)
    with open("average_mean_high_skew_output.json", "w") as f:
        json.dump(skewed_but_average,f)

    return(enhancements, CME_or_lag, swab_high_mean, skewed_but_average)








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
            elif single_mean < 1.4:
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
    #print(len(means))

    with open("all_means_output.json", "w") as f:
        json.dump(means, f)

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


