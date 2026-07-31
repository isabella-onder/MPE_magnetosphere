from astropy.io import fits
import numpy as np
import scipy
import scipy.stats
import matplotlib.pyplot as plt
import json
import glob
import ast



#extract all the files
with open('all_flared_obsIDs.json','r') as f:
    #obsID_enhancement_list = json.load(file)
    obsID_enhancement_list = [line.strip() for line in f if line.strip()]
    #print(obsID_enhancement_list, type(obsID_enhancement_list))

#extract the big list with ALL the info (ObsID with corresponding timestamps (timestamps of flares detected in eROSITA i.e. passed the bayesian blocks threshold but but in xmm time) and rev): indices and order corespond between Obs_ID list and this list with all
with open('big_list.json','r') as file:
    big_list = file.read()
    big_list = ast.literal_eval(big_list)
    #print(big_list, type(big_list))

#extracting ALL the event count rates for eROSITA and corresponding times
with open('all_times.json', 'r') as f:
    all_erosita_times = f.read()
    all_erosita_times = ast.literal_eval(all_erosita_times)
with open('all_normalised_n.json', 'r') as f:
    all_normalised_n = f.read()
    all_normalised_n = ast.literal_eval(all_normalised_n)

def position_extracter(filepath,starting_time, ending_time):
    #print(filepath)

    #need to adapt the filepath for background flare such that it is indeed the orbit file
    filepath_orb = filepath.rsplit('/', 2)[0]
    filepath_to_search = filepath_orb + "/PN*/*ORB*.FIT*"
    #print(filepath_to_search, 'filepath_to_search')
    filepath_orb = glob.glob(filepath_to_search)[0]

    #print('this is the new path file', filepath_orb, type(filepath_orb))



    with fits.open("/data36s/bella/erodat/erosita/omnia_orbit_merged.fits") as hdul_ero, fits.open(filepath_orb) as hdul_xmm:
    
        #we look at positions at both starting and ending time of each obs_ID (erosita will stary roughly the same, simply to also track xmm
        

        # coords for erosita

        data_erosita = hdul_ero[1].data
        eroday = data_erosita['eROday']
        time_erosita = data_erosita['TIME']
        x_ecl_erosita = data_erosita['X_ECL']
        y_ecl_erosita = data_erosita['Y_ECL']
        z_ecl_erosita = data_erosita['Z_ECL']
        #print(len(x_ecl_erosita), 'this is len x_ecl_erosita')
        #need to convert the starting time back to erosita timestamps to be able to compare
        starting_time_ero = starting_time - 63061262 #this is the same time delay as used in comparator
        ending_time_ero = ending_time - 63061262

        #index of starting time row
        idx_temp = np.searchsorted(time_erosita, starting_time_ero)
        idx = idx_temp - 1 #since searchsorted gives where it would be placed, go to previous row for previous stamp
        
        #print('this is the eroday at that index - 1', eroday[idx-1])
        #print('this is the starting time of the eroday', starting_time_ero, 'and the corresponding time_erosita close', time_erosita[idx-1])
        #print('this is idx from search sorted',idx)
        x_ecl_ero_start = x_ecl_erosita[idx]
        y_ecl_ero_start = y_ecl_erosita[idx]
        z_ecl_ero_start = z_ecl_erosita[idx]

        idx_end_temp = np.searchsorted(time_erosita, ending_time_ero)
        idx_end = idx_end_temp - 1
        x_ecl_ero_end = x_ecl_erosita[idx_end]
        y_ecl_ero_end = y_ecl_erosita[idx_end]
        z_ecl_ero_end = z_ecl_erosita[idx_end]

        #print('just to check idxs', idx, idx_end)
        erosita_coords = [(x_ecl_ero_start, y_ecl_ero_start, z_ecl_ero_start), (x_ecl_ero_end, y_ecl_ero_end, z_ecl_ero_end)]



        # coords for xmm
        data_xmm = hdul_xmm[1].data

        time_xmm = data_xmm['TIME']
        x_gei_xmm = data_xmm['GEI_X']
        y_gei_xmm = data_xmm['GEI_Y']
        z_gei_xmm = data_xmm['GEI_Z']

        x_gse_xmm = data_xmm['GSE_X']
        y_gse_xmm = data_xmm['GSE_Y']
        z_gse_xmm = data_xmm['GSE_Z']

        #index of starting/ending rows
        idx_xmm = np.searchsorted(time_xmm, starting_time)
        
        x_gei_xmm_start = x_gei_xmm[idx_xmm]
        y_gei_xmm_start = y_gei_xmm[idx_xmm]
        z_gei_xmm_start = z_gei_xmm[idx_xmm]

        x_gse_xmm_start = x_gse_xmm[idx_xmm]
        y_gse_xmm_start = y_gse_xmm[idx_xmm]
        z_gse_xmm_start = z_gse_xmm[idx_xmm]

        idx_end_xmm = np.searchsorted(time_xmm, ending_time)
        
        x_gei_xmm_end = x_gei_xmm[idx_end_xmm]
        y_gei_xmm_end = y_gei_xmm[idx_end_xmm]
        z_gei_xmm_end = z_gei_xmm[idx_end_xmm]

        x_gse_xmm_end = x_gse_xmm[idx_end_xmm]
        y_gse_xmm_end = y_gse_xmm[idx_end_xmm]
        z_gse_xmm_end = z_gse_xmm[idx_end_xmm]

        #print('idxs for xmm', idx_xmm, idx_end_xmm)
        xmm_gei_coords = [(x_gei_xmm_start, y_gei_xmm_start, z_gei_xmm_start), (x_gei_xmm_end, y_gei_xmm_end, z_gei_xmm_end)]
        xmm_gse_coords = [(x_gse_xmm_start, y_gse_xmm_start, z_gse_xmm_start), (x_gse_xmm_end, y_gse_xmm_end, z_gse_xmm_end)]

        #converting xmm gei to ecliptic such as to be able to compare with erosita
        xmm_ecl_coords = []
        eta = np.deg2rad(23.5)
        for coords in xmm_gei_coords:
            v = np.array(coords)
            M = np.array([
                [1, 0, 0],
                [0, np.cos(eta), np.sin(eta)],
                [0, - np.sin(eta), np.cos(eta)]
            ])
            ecl = M @ v
            xmm_ecl_coords.append(ecl)
        #print(xmm_gei_coords, xmm_ecl_coords)
    return erosita_coords, xmm_gei_coords, xmm_gse_coords, xmm_ecl_coords









#lightcurve_plotter: cycles through all the filepaths in obsID_enhancement_list (from comparator.py most likely) and produces plot. For now, input manually
def lightcurve_plotter(obsID_enhancement_list, big_list, all_erosita_times, all_normalised_n, plot_validates):
    np.set_printoptions(threshold=np.inf)
    #print('this is obsID_enhancement_list at the cusp of the for loop in lightcurve_plotter', obsID_enhancement_list)
    #print('and this is the corresponding obs_entry info from big_list', big_list)
    correlator_big_list = []
    z_on_same_side = []
    for filepath, obs_entry in zip(obsID_enhancement_list, big_list):
        #print('this is obs_entry', obs_entry)
        correlator_small_list = []
        correlator_small_list.append(obs_entry[0]) #adding the obsID number to the small list
        title = obs_entry[0]
        if filepath == 'file_did_not_exist':
            print('this filepath did not exist')
            continue
        print('this filepath is about to be plotted', filepath)
        with fits.open(filepath) as hdul: #opening the xmm obs_ID files to extract the already processed background rates and corresponding timestamps
            time = hdul[1].data['TIME']
            time_start = min(time)
            time_end = max(time)
            #print('we need to check the ranges, especially time start: might have to put an error pass in case it s null', time_start, time_end)
            print('this is a time step for xmm:', time[1]-time[0])
            

            #this was the name to extract rates from the original files with different count rates
            #rate = hdul[1].data['RATE']
            #this is from the files with bins of 120s recreated by MJF
            rate = hdul[1].data['COUNTS']

            timestamps = obs_entry[1]
            lengths = obs_entry[2]
            
            indices_erosita = [i for i, t in enumerate(all_erosita_times) if time_start < t < time_end] #might have to make them all integers
            
            #I want it to take a few extra minutes off at the end rather than right at the start (since more interesting to shift that way: so add say 5 so that it decales by 120 * 5 seconds (takes 5 extra readings from after the end of the ObsID for erosita
            #indices_erosita = [i + 5 for i in indices_erosita]

            #print(indices_erosita, len(indices_erosita))
            #print('len(all_normalised_n)', len(all_normalised_n), 'len(all_erosita_times)', len(all_erosita_times))
            appropriate_erosita_rates = [all_normalised_n[i] for i in indices_erosita]
            appropriate_erosita_times = [all_erosita_times[i] for i in indices_erosita]
            #print('this is the start and end of appropriate erosita times', appropriate_erosita_times[0], appropriate_erosita_times[-1])
            #print('and these are the rates during those times', appropriate_erosita_rates)
            #plt.plot(time, rate, marker = '.')

            #for now is in a 'try' loop in case there is error due to there not being such files in MOS
            try:
                erosita_coords, xmm_gei_coords, xmm_gse_coords, xmm_ecl_coords = position_extracter(filepath, time_start, time_end)
                
                #want to check whether during an obs_ID, SRG and XMM are on the same side of the ecliptic or not
                #if erosita_coords[0][3] * xmm_ecl_coords[3] > 1: #if they are, will be the same side of z_ecl
                #    z_on_same_side.append(True)
                #elif erosita_coords[0][3] * xmm_ecl_coords[3] <1: #if not, will not be
                #    z_on_same_side.append(False)
                #else: #i.e. necessarily 0
                #    z_on_same_side.append(np.nan)
                #print(i'are these some appropriate erosita_coords?', erosita_coords, xmm_gei_coords, xmm_gse_coords)

            except:
                #z_on_same_side.append(np.nan)
                print('there was no orb file: a nan was added to the z comparator')

            #unpacking the coordinates to print
            xmm_ecl_start_x = xmm_ecl_coords[0][0]
            xmm_ecl_end_x = xmm_ecl_coords[1][0]
            xmm_ecl_start_y = xmm_ecl_coords[0][1]
            xmm_ecl_end_y = xmm_ecl_coords[1][1]
            xmm_ecl_start_z = xmm_ecl_coords[0][2]
            xmm_ecl_end_z = xmm_ecl_coords[1][2]

            xmm_ecl_start = list(xmm_ecl_coords[0])
            xmm_ecl_end = list(xmm_ecl_coords[1])
            erosita_ecl_start = list(erosita_coords[0])    

            
            #want to check whether during an obs_ID, SRG and XMM are on the same side of the ecliptic or not
            if erosita_ecl_start[2] * xmm_ecl_start_z > 1: #if they are, will be the same side of z_ecl
                z_on_same_side.append([True, title])
            elif erosita_ecl_start[2] * xmm_ecl_start_z <1: #if not, will not be
                z_on_same_side.append([False, title])
            else: #i.e. necessarily 0
                z_on_same_side.append([np.nan, title])

            #plotting(perhaps)
            if plot_validates:
                fig, ax1 = plt.subplots()

                ax1.set_xlabel("Time (xmm seconds)")
                ax1.set_ylabel("Counts/sec XMM")
                ax1.plot(time, rate, marker = '.', label = 'xmm', color = 'royalblue', ms = 2)

                ax2 = ax1.twinx()
                ax2.set_ylabel("Counts/sec eROSITA")
                ax2.plot(appropriate_erosita_times, appropriate_erosita_rates, marker = '.', ms = 2, color = 'chocolate', label = 'erosita')

            #this puts some background colour on the bayesian blocks that were deemed to be flares i.e. flagged since above a threshold
            #green if long, pink if shorter than 100s
                for start, duration in zip(timestamps, lengths):
                    if duration > 100:
                        color = 'chocolate'
                        alpha = 0.3
                    else:
                        color = 'wheat'
                        alpha = 0.9

                    end = start + duration
                    #print('start', start, 'end', end, 'duration', duration)
                    plt.axvspan(start, end, color = color, alpha = alpha)

                txt = 'Satellite eclipic coords during this obsID'
                txt_xmm = 'XMM coords start:' + str(xmm_ecl_start) + '\n XMM coords end: '+str(xmm_ecl_end)
                txt_ero = 'eROSITA start: '+ str(erosita_ecl_start)

                fig.text(0.5, 0.04, txt, ha='center')
                fig.text(0.5, 0.02, txt_xmm, ha='center', color = 'royalblue')
                fig.text(0.5, 0.01, txt_ero, ha='center', color = 'chocolate')

                plt.xlabel("Time (xmm seconds)")
                plt.title('Plotted obsID ' + title)
                fig.legend()
                plt.show()

            #print('xmm times', time)
            #print('xmm rates', rate)
            #print('erosita times', appropriate_erosita_times)
            #print('erosita rates', appropriate_erosita_rates)
            
            #adding all rates and timestamps to the small correlator list (i.e. all the info for THIS obsID)
            correlator_small_list.append(time.tolist())
            correlator_small_list.append(rate.tolist())
            correlator_small_list.append(appropriate_erosita_times)
            correlator_small_list.append(appropriate_erosita_rates)

            #adding the small list to the large, overarching one
            correlator_big_list.append(correlator_small_list)
            
            #unpacking the coordinates to print
            xmm_ecl_start_x = xmm_ecl_coords[0][0]
            xmm_ecl_end_x = xmm_ecl_coords[1][0]
            xmm_ecl_start_y = xmm_ecl_coords[0][1]
            xmm_ecl_end_y = xmm_ecl_coords[1][1]
            xmm_ecl_start_z = xmm_ecl_coords[0][2]
            xmm_ecl_end_z = xmm_ecl_coords[1][2]

            xmm_ecl_start = list(xmm_ecl_coords[0])
            xmm_ecl_end = list(xmm_ecl_coords[1])
            erosita_ecl_start = list(xmm_ecl_coords[0])

    
        #return time, rate, appropriate_erosita_times, appropriate_erosita_rates
    with open("correlator_FINAL.json", "w") as f:
        json.dump(correlator_big_list, f, indent = 2)
        print('the correlator_big_list containing ObsID, timestamps and rateswas succesfully dumped')
        print('z_on_same_side test', z_on_same_side)

lightcurve_plotter(obsID_enhancement_list, big_list, all_erosita_times, all_normalised_n, True)
