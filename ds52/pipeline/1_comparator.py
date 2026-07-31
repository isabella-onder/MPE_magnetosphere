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
from scipy.ndimage import label
from collections import defaultdict

import glob

plot_validity = False

with open('all_enhancements_output.json','r') as file:      #loading all the eROdays which have enhancements, as per detected by machine.py (chi_squared.py)
	eRO_enhancements_list = json.load(file)
#print(eRO_enhancements_list)



def lightcurve_plotter(eROday, energy_threshold, binsize): #usually we will be inputting eRO_enhanements_list, but will have the flexibility to change threshold and binning size
    filename = "e_6_" + str(eROday) + "_002_c030.fits"
    filepath = "/data36s/bella/erodat/erosita/eROday/" + filename                                                                                                                                                                                                                                                                                                                                                                

    #plot_validity = True
    with fits.open(filepath) as hdul_ero, fits.open("/data36s/bella/erodat/erosita/artxc_catalogue_ecliptic.fits") as hdul_art, fits.open("/data36s/bella/erodat/erosita/maxi_bright_src.fits") as hdul_maxi:
        print('this eROday file is indeed being read', filename)


        non_filtered_events = hdul_ero[1].data                                                                          
        events = non_filtered_events[non_filtered_events['PI'] > energy_threshold]   
        #print('this is the len of events', len(events))
	
	#just in case there are no registered events then terminating code
        if len(events)==0:
            print('Len(events)==0. Something has gone very wrong with eROday'+str(eROday))
            mean, median, mode, skew = -1000, -1000, -1000, -1000 #putting 1000 as an absurd value to flag it and leave it
            return mean, median, mode, skew, -1000, -1000


	# applying the mask
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
	
	#extracting all counts that are not within the mask setting the time to the scale of an eROday (to make the maths easier - can then reference "time" when need the actual time)
        time = filtered_events['TIME']
        reset_time = time - min(time)
        bin_time = binsize
        poubelles = int((max(time)-min(time))/bin_time)
        n,bins = np.histogram(reset_time, poubelles, (min(reset_time),max(reset_time)))
        n, bins_non_reset = np.histogram(time, poubelles, (min(time), max(time)))
        normalised_n = n/bin_time
        time_in_hours = bins/3600

        bin_nb_for_histogram = 12

        #removing all zeroes to avoid changing the stats
	#NB: even though the zeroes still appear on the stairs plot, they are NOT included in the histogram and statistical parameters' calculationsi
        idx_zeroes = np.where(normalised_n == 0)[0]
        idx_zeroes_and_around = list(idx_zeroes)
        for index in idx_zeroes:
            idx_zeroes_and_around.append(index-1)
            idx_zeroes_and_around.append(index+1)
        idx_zeroes_and_around = list(set(idx_zeroes_and_around))

        normalised_n_without_zeroes = np.delete(normalised_n,idx_zeroes_and_around)
        #reset_time_without_zeroes = np.delete(reset_time, idx_zeroes_and_around)
        #print('checking the lengths of reset time', len(reset_time), 'and without zeroes', len(reset_time_without_zeroes))
        time_in_hours_without_zeroes = np.delete(time_in_hours, idx_zeroes_and_around) ### check whether I detele zero times or whether I just use counts that do not have zeroes removed
        #print(time_in_hours, time_in_hours_without_zeroes, 'the times vs times without zeroes and lengths respectively', len(time_in_hours), len(time_in_hours_without_zeroes))
       
       #just in case there were only 0 then terminating code
        if len(normalised_n_without_zeroes)==0:
            print('Len(normalised_events_without_zeroes)==0. There are NO counts once stars are substracted - that is suspicious, check len(events) - eROday'+str(eROday))
            mean, median, mode, skew = -1000, -1000, -1000, -1000 #putting 1000 as an absurd value to flag it and leave it
            return mean, median, mode, skew, -1000, -1000
	   
        
       #for plotting, when necessary (removed since do not want it to appear every iteration
        if plot_validity == True:
            figue, axs = plt.subplots(nrows = 1, ncols = 2)
            axs[0].stairs(normalised_n_without_zeroes, time_in_hours_without_zeroes, color = 'chocolate')
            axs[1].stairs(normalised_n, time_in_hours, color = 'chocolate')
            figue.supylabel('Count rate')
            figue.supxlabel('Time in hrs')
            figue.suptitle('eROSITA detections \n Raw binning on the left, zeroes removed on the left')
            #plt.show()

        #parameters WITHOUT the zeroes
        mean= np.mean(normalised_n_without_zeroes)
        median = np.median(normalised_n_without_zeroes)
        mode = scipy.stats.mode(normalised_n_without_zeroes)
        skew = scipy.stats.skew(normalised_n_without_zeroes)
        return(mean,median,mode,skew,normalised_n, normalised_n_without_zeroes, time, time_in_hours_without_zeroes, bins_non_reset)
        print('Mean: ',mean,'\nMedian: ',  median,'\nMode: ', mode[0][0],'\nSkew: ',skew)  
	
def flare_finder(eROday_list, energy_threshold, bins):
    flare_ranges = []
    depressions = []
    all_xmm_time_edges = []
    xmm_times_flares_start_all = []
    xmm_times_flares_lengths_all = []
    all_normalised_n = []
    all_times = []

    #plot_validity = True

    #different threshold for different bin sizes, to account for dilution
    if bins == 32:
        skew_threshold = 0.4
        #print('the bins are 32, the chosen skew threshold was 0.4')

    elif bins == 120:
        skew_threshold = 0.7
        #print('the bins are 120, the chosen skew threshold was 0.7')
	
	
    else: #default
        skew_threshold = 0.7
        #print('the default skew threshold was chosen')

    

    for eROday in eROday_list:
        #print('this is the eROday', eROday)

        #accomodating for the different setups and therefore different background rates
        if (eROday > 43411 and eROday < 46309):         
            mean_threshold = 1.801

        else:
            mean_threshold = 2.001

    
    #retrieving the data for every eROday: stats, counts, corresponding times. Setting parameters for background vs enhanceement thresholds accordingly
        mean, median, mode, skew, normalised_n_full, normalised_n, time, time_in_hours, erosita_bins = lightcurve_plotter(eROday, energy_threshold, bins)
        #print('this is the output of lightcurve_plotter ', mean, median, mode, skew, normalised_n, time_in_hours)
        bg = median
        std = np.std(normalised_n)
        #print('this is bg, std, threshold: ', bg, std, threshold)


        #for high skew high mean: finds by looking at 3std deviation from bg median to annonate peaks (finds end and start bin and adds to list)
        if (mean > mean_threshold and skew > skew_threshold):
            threshold = bg + 1.3 * std
            #this is the strong type: both high mean and skew
            #need only one std deviation - the peaks are so massive anyways and std large that anything deviating from the median probably is a peak

            print('the skew was high the mean too')
            
            flare_mask = normalised_n > threshold
            labeled, n_flares = label(flare_mask)
            times_to_plot = []
            times_to_plot_end = []
            for i in range(1, n_flares+1):
                idx = np.where(labeled == i)[0]
                start, end = time[idx[0]], time[idx[-1]]   #ask Michael and just add on the correct number of seconds
                times_to_plot.append(time_in_hours[idx[0]])
                times_to_plot_end.append(time_in_hours[idx[-1]])
                flare_ranges.append([start, end])
            
        #for high mean, low skew (these went through swab)
        elif (skew < skew_threshold and mean > mean_threshold): 
            print('the skew was low, the mean high')
            threshold = bg + 1 * std            
            #these are the ones that went in for a swab: low skew, high mean
            if skew < 0:
                threshold = 1000 #if there is a negative skew, it is a depression rather than a feature - will save to the side but do not want it to spoil data
                depressions.append(eROday)


        elif (skew > skew_threshold and mean < mean_threshold): #these were the potential sources or at least peaks: high skew, low mean
            print('the skew was high the mean low')
            threshold = bg + 3 * std
            flare_mask = normalised_n > threshold
            labeled, n_flares = label(flare_mask)
            times_to_plot = []
            times_to_plot_end = []
            for i in range(1, n_flares+1):
                idx = np.where(labeled == i)[0]
                start, end = time[idx[0]], time[idx[-1]]   #ask Michael and just add on the correct number of seconds
                times_to_plot.append(time_in_hours[idx[0]])
                times_to_plot_end.append(time_in_hours[idx[-1]])
                flare_ranges.append([start, end])



        else:
            threshold = bg + 1 * std

        
        #print('this is bg, std, threshold', bg, std, threshold)
        


        time_reset = time - min(time)
        #print('this is min(time)', min(time), 'this is min(time_reset)', min(time_reset))
        edges = bayesian_blocks(time_reset, fitness='events', gamma= 0.000001) #letting it do the binning
        #edges = edges - edges[0]
    
        
        if plot_validity:
            #plot to see where the delimitations are
            fig_bb, ax = plt.subplots()
            ax.stairs(normalised_n, time_in_hours, color = 'chocolate')
            for e in edges:
                ax.axvline(e/3600, color='black', ls='--')
            fig_bb.supxlabel('Time in hours')
            fig_bb.supylabel('Count rate')
            fig_bb.suptitle('eROSITA lightcurve for eROday' +str(eROday)+'\n Vertical lines depict bayesian block splitting')
            #fig_bb.show()
        
       #extracting the number of counts in each bayesian block
        new_n, _ = np.histogram(time_reset, edges)
        
        #unncessary: was comparing with astropy.stats which just uses wrapper around bayesian_blocks and cannot modify precision    
        '''astro_new_n, astro_bins_edges = astropy.stats.histogram(time_reset, bins = 'blocks')
        plt.stairs(normalised_n, time_in_hours)
        for a in astro_bins_edges:
            plt.axvline(a/3600, color='blue', ls='--')
        plt.show()
        '''
    
        #finding the count rates in each block 
        diffs_edges = [edges[i] - edges[i-1] for i in range(1, len(edges))]
        normalised_new_n = [n / d for n,d in zip(new_n, diffs_edges)] #these should be the average count/sec in every block
        #print('Edges: ', edges, len(edges), '\n New_n (nb of events in each block): ', new_n, len(new_n), '\n Length of each block: ', diffs_edges, len(diffs_edges), '\n Normalised count rate (ct/sec): ', normalised_new_n, len(normalised_new_n))
        #looping through the bayesian block rates to remove ghost spikes: i.e. extremely narrow, less than a second dozen counts which threfore produce very high rates and outweigh the rest
        ghost_spikes = []
        for i in range(len(normalised_new_n)):
            if normalised_new_n[i] > 15 and diffs_edges[i] < 1.0:
                ghost_spikes.append((edges[i] + min(time), normalised_new_n[i]))
                normalised_new_n[i] = 0
                print('a ghost spike has been detected:', ghost_spikes)

        if plot_validity:
            fig_BB, axs_BB = plt.subplots()
            #visually showing the bayesian blocks, grouped: count rate per block over time
            axs_BB.stairs(normalised_new_n, edges, color = 'chocolate')
            fig_BB.suptitle('Bayesian blocks stairs: counts are grouped by block and depict average count rate of the block')
            fig_BB.supxlabel('Time in hours')
            fig_BB.supylabel('Average count rate in BB')
            #plt.show()

        #making note of the indices where normalised_new_n was above our mean threshould for flares
        arr = np.array(normalised_new_n) #for np.where to work it needs an array it seems
        high_rate_bins_indices = np.where(arr > mean_threshold)[0] #returns array of arrays in case there are multiple variables/conditions being checked



        #converting eROday timestamps to XMM time
        indices = [np.argmin(np.abs(time_reset - e)) for e in edges]
        #print('this is the indices info', indices, len(indices))
        extract_time = [int(time[i]) for i in indices] #no need for accuracy to less than unit
        #print(extract_time, 'this is extract time')
        xmm_time = [ x + 63061262 for x in extract_time] #there may be shift up to 1.3 seconds
        
        xmm_times_flare_start = [xmm_time[i] for i in high_rate_bins_indices] #getting the starting time of all flares (i.e. above the mean threshold)
        xmm_times_flares_start_all = xmm_times_flares_start_all + xmm_times_flare_start
        corresponding_lengths = [diffs_edges[i] for i in high_rate_bins_indices] #getting the corresponding length in seconds of each flare
        xmm_times_flares_lengths_all = xmm_times_flares_lengths_all + corresponding_lengths
        all_xmm_time_edges = xmm_time + all_xmm_time_edges #adding it to an overarching, larger list
        #print('did depression hit?', depressions)
    

        #extracting the actual counted event rates, normalised_n, and converting all corresponding times from eROSITA to XMM to be able to plot them together
        time = list(time)
        erosita_bins = list(erosita_bins)
        erosita_bins.pop()
        #normalised_n = list(normalised_n)
        normalised_n = list(normalised_n_full)
        #print('these are the lengths before any conversion', len(normalised_n), len(erosita_bins))
        all_normalised_n = all_normalised_n + normalised_n
        xmm_times = [t + 63061262 for t in erosita_bins] #converting all the times
        all_times = all_times + xmm_times #adding to a large, overarching list
        #print('these are the lengths of the overarching lists', len(all_normalised_n), len(all_times))
        
        


        
        ''' #greyed it out since bayesian blocks is more precise
        #'try' so that it does it if and only if peaks were detected (low mean high skew) i.e. if times_to_plot exists
        try:
            arbitrary = [1] * len(times_to_plot)
            plt.scatter(times_to_plot, arbitrary, s = 20,color = 'r')
            plt.scatter(times_to_plot_end, arbitrary,s = 15, color = 'blue')
            plt.stairs(normalised_n, time_in_hours)
            print('times to plot are', times_to_plot)

            plt.title('High skew and low mean eROday - detected peaks')
            plt.show()

        except UnboundLocalError:
    
        print('UNBOUND LOCAL ERROR')
            continue
         '''   

        #print('the returned xmm_times_flare_start', xmm_times_flare_start)
        plt.show()
    #print('final all_xmm_time_edges', all_xmm_time_edges, len(all_xmm_time_edges))
    #print('final xmm_times_flares_start_all',xmm_times_flares_start_all, len(xmm_times_flares_start_all))
    #print('final xmm_times_flares_lengths_all', xmm_times_flares_lengths_all, len(xmm_times_flares_lengths_all))
    #print('these are the final lengths of the overarching erosita data lists', len(all_normalised_n), len(all_times))
    
    #print('these are the supposedly problematic lengths: len(all_normalised_n)', len(all_normalised_n), ' len(all_times)', len(all_times))

    #writing all the counts and corresponding times studied above
    with open("all_normalised_n.json", "w") as f:
        f.write(str(all_normalised_n))
    with open("all_times.json", "w") as f:
        f.write(str(all_times))


    return xmm_times_flares_start_all, xmm_times_flares_lengths_all #returns the start and duration of the flares detected in erosita, but has already been translated to xmm time
	                   
def xmm_flare_finder(eROday_list, threshold, bins): # can put it in the enhancement eROday list eventually
    xmm_timestamps, xmm_corresponding_lengths = flare_finder(eROday_list, threshold, bins)
    with fits.open("/data36s/bella/erodat/xmm/xsa_archive_2019_2022_repaired.fits") as hdul:
        data = hdul[1].data

        obs_ID = data['OBSERVATION']
        xmm_times = data['TIME_0_XMM']
        revs = data['REVOLUTION']
        xmm_end_times = data['TIME_E_XMM']

        xmm_times_sorted = np.array(xmm_times) #sorted for searchsorted to work
        
        #extracting the indices of the obsIDs of the timestamps we study, to then be able to look at rows more easily - had to change method to avoid timestamp being put in the adjacent obs_id, since edge is closer
        rows = np.searchsorted(xmm_times, xmm_timestamps, side='left')
        #print(rows, 'and in comparison indices_test')
        indices = [row - 1 for row in rows]
        obs_IDs_repeat = []
        for i in indices:
            if obs_ID[i][0] == 9:
                obs_IDs_repeat.append(obs_ID[i])
            else:
                obs_IDs_repeat.append(obs_ID[i])

        #print('this is the list of obs_IDs repeat, hopefully corresponding to every event', obs_IDs_repeat)

        #the starting times of each ObsID
        corresponding_times = [xmm_times[i] for i in indices]
        corresponding_times_set = list(set(corresponding_times))
        corresponding_times_set.sort()
        #print('these are the corresponding_times (i.e. what should be the starting time of each obs_ID', corresponding_times)
        corresponding_end_times = [xmm_end_times[i] for i in indices]
        corresponding_end_times_set = list(set(corresponding_end_times))
        corresponding_end_times_set.sort()
        #print('these are the corresponding_end_times (i.e. what should be the ending point of each obs_ID', corresponding_end_times)

        corresponding_revs_repeat = [revs[i] for i in indices]

        #creating big list of format list =  [obs_ID_1, [timestamps_1_1, timestamp_1_2, ...]], [obs_ID_2, [timestamp_2_1, timestamp_2_2, ...]], ...]
        #here, timestamp refers to the starting timestamp of a bayesian block which was flagged as a flare, with respective durations
        groups = defaultdict(lambda: {'timestamps': [], 'revs': set(), 'lengths': []})
        for t, obs, rev, length, end in zip(xmm_timestamps, obs_IDs_repeat, corresponding_revs_repeat, xmm_corresponding_lengths, corresponding_end_times):
            if t > end:
                #print('the t > end was initiated, hopefully it has continued for t ', t, ' obs ' , obs)
                continue
            else:
                groups[obs]['timestamps'].append(t)
                groups[obs]['revs'].add(rev)
                groups[obs]['lengths'].append(length)
        complete_obsID_timestamp = [[obs, data['timestamps'],data['lengths'], list(data['revs'])[0]] for obs, data in groups.items()]
        #print('my big list: obs ID with corresponding timestamps and rev',complete_obsID_timestamp)
        

        #relevant_obs_ID = list(set([int(obs_ID[i]) for i in indices]))
        corresponding_rev = list(set([int(revs[i]) for i in indices]))
        corresponding_rev.sort()

        #creating similar big list where the timestamps are relative to the obs_ID (i.e. t = 0 at the start of each obsID)
        complete_obsID_timestamps_normalised = []
        for obs_entry, starting_time in zip(complete_obsID_timestamp, corresponding_times_set):
            obs_id = obs_entry[0]
            rev = obs_entry[-1]
            #print('this is the obs_entry', obs_entry[0], 'this is the obs_id', obs_id, 'this is the corresponding starting_time', starting_time)
            new_timestamps = [t - starting_time for t in obs_entry[1]]
            complete_obsID_timestamps_normalised.append([obs_id, new_timestamps, rev])
        #print('my other big list: obs ID with corresponding timestamps in seconds starting from obs ID start', complete_obsID_timestamps_normalised)

        
        '''
        str_revs = [str(day) for day in corresponding_rev]
        str_obs_ID = [str(0) + str(obs) for obs in relevant_obs_ID]

        with open("list.txt", "a") as f:
            for day in str_revs:
                for obs_ID in str_obs_ID:
                    print(day, obs_ID, 'the file looked for this day and obs_ID')
                    files = glob.glob(f"/xmm/archive/products/{day}/{obs_ID}/PN*/P*PN*FBKTSR*.FIT*")
                    print(files)
                for filename in files:
                    f.write(filename + "\n")
        '''

        with open("all_flared_obsIDs.json", "w") as f:
            data = []
            for obs_entry in complete_obsID_timestamps_normalised:
                obs_id = obs_entry[0]
                #obs_id_slew = '9' + obs_id[1:]
                #print('this is obs_ID_slew', obs_id_slew)
                rev = obs_entry[2]
                files = glob.glob(f"/xmm/archive/products/{rev}/{obs_id}/PN*/P*PN*FBKTSR*120*.FIT*")
                print('the file looked for this rev and obs_ID', rev, obs_id)
                #if files == []:
                #    continue
                #data.append(files)

                #if not files: #this obs_id may have been a slew observation isntead
                    #files = glob.glob(f"/xmm/archive/products/{rev}/{obs_id_slew}/PN*/P*PN*FBKTSR*.FIT*")
                    #print('a slew observation was searched for instead')
                

                if not files: #i.e. if there was no PN probably - make it look in MOS1 folder instead
                    files = glob.glob(f"/xmm/archive/products/{rev}/{obs_id}/M1*/P*M1*FBKTSR*120*.FIT*")
                    #print('the mos file was searched for')
                    #print('and here is the files variable', files)

                #if not files: #now check MOS for a slew observation
                    #files = glob.glob(f"/xmm/archive/products/{rev}/{obs_id_slew}/M1*/P*M1*FBKTSR*.FIT*")

                if not files: 
                    #print('the if not files statement has been iterated twice (there was neither PN nor MOS1)')
                    f.write("file_did_not_exist \n")
                    continue

            #print('this is hopefully all files_list', data)
                for filename in files:
                    print('this file was added', filename)
                    f.write(filename + "\n")
                #data.append(filename)
                #json.dump(data,f, indent=4)

        with open("big_list.json", "w") as f:
            f.write(str(complete_obsID_timestamp))
            #json.dump(complete_obsID_timestamp, f, indent=4)

    
           
