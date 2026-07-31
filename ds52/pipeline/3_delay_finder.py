import numpy as np
import usefulness as u
import matplotlib.pyplot as plt
import json

plot_wanted = False

#extracting the info output from xmm_plotter: a large list where each element is a list containing ['obsID number', [xmm times], [xmm rates], [erosita times], [erosita rates]]
with open("correlator.json", "r") as f:
    large_list = json.load(f)


#I should be able to have all the functions outside of my loop hopefully

def xcorr_data(d1, d2, maxlags=None, normed=True):
    # IS united tested: returns the same thing as plt.xcorr!
    # Alternative code lifted from LIGO analysis
    # My version is cross_correlate
    d1 = np.asarray(d1)
    d2 = np.asarray(d2)

    #print('these are the lengths of d1 and d2 respectively in xcorr_data', len(d1), len(d2))
    N = len(d1)

    if maxlags is None:
        maxlags = N - 1

    lags = np.arange(-maxlags, maxlags + 1)
    c = np.correlate(d1 - np.mean(d1), d2 - np.mean(d2), mode='full')

    if normed:
        scale = np.std(d1) * np.std(d2) * N
        c = c / scale

    mid = len(c) // 2
    start = mid - maxlags
    stop = mid + maxlags + 1
    delays = lags / f
    return delays, c[start:stop]


def cross_correlate(d_1:np.array, d_2:np.array, t_1:np.array, t_2:np.array, use_mathematically_correct_normalization=False):
    if not np.all(t_1 == t_2):
        raise ValueError("The two time arrays must be the same.")
    t = t_1



    d_1 = d_1 - np.mean(d_1)
    d_2 = d_2 - np.mean(d_2)

    #print(i'these are the lengths of d1 and d2 respectively in cross_correlate', len(d_1), len(d_2))
    assert len(d_1) == len(d_2), "The two data arrays must have the same length."

    n = len(d_1)
    #l = np.linspace(-(n-1), (n-1), 2*(n-1)+1)
    l = np.linspace(0, (n-1), (n-1))        #look only at positive delays, since anything else is unphysical (occurs first in XMM)
    f = 1/(t[1]-t[0])
    taus = l/f

    # loop for all l's:

    correlations = []
    for l_star in l:
        l_star = int(l_star)
        # loop for a single l_star:
        summands = []
        for i in range(int(n-np.abs(l_star))):  # upper limit: n-l-1
            summands.append(d_1[i] * d_2[i+l_star])

        if use_mathematically_correct_normalization:
            normalization = n-np.abs(l_star)  # mathematically correct
        else:
            normalization = n  # suppresses regions of large delays, which we are not interested in because there cant be physically sourced correlations
            # after dividing with the noise power spectrum estimate that are on the order of ~10ms! These would be random noise correlations ... I think
        res_for_l_star = sum(summands) / normalization
        correlations.append(res_for_l_star)

    correlations = np.array(correlations)/(np.std(d_1)*np.std(d_2))   # normalize with the standard deviations if you want
    return taus, correlations  # return time delays and C(τ)

def divisor(t: np.array, x: np.array, window_length: float):
    """
    Subdivides a given (x, t) time series into non-overlapping windows of duration `window_length` (in seconds).
    Assumes t is sorted and uniformly sampled.

    :param x:               The time series values.
    :param t:               The time stamps (1D array, same length as x).
    :param window_length:   Window length in seconds.
    :return:                Tuple of (windows of t, windows of x), both shaped (num_windows, samples_per_window)
    """
    x = np.array(x)
    t = np.array(t)

    dt = t[1] - t[0]  # assume uniform sampling
    samples_per_window = int(window_length / dt)

    total_samples = len(t)
    usable_samples = (total_samples // samples_per_window) * samples_per_window

    x = x[:usable_samples]
    t = t[:usable_samples]

    res_x = x.reshape(-1, samples_per_window)  # this must have been chat gpt I don't use reshpae
    res_t = t.reshape(-1, samples_per_window)

    return res_t, res_x

overarching_obsIDs = []
more_than_2stds = []
more_than_3stds = []
overarching_optimal_delays = []
#looping through each element of the obsIDs, to correlate each
for i in range(len(large_list)):
    obsID_data = large_list[i] #just using the first obsID as the timestamps match
    obsID, time_xmm, xmm, time_erosita, erosita = obsID_data #unpacking the list (each unpacked list was an element of the large overarching list)
    try:
        start = max(time_xmm[0], time_erosita[0])
        end = time_erosita[-1]
    except IndexError:
        print('there was an index error for time_xmm[0] and time_erosita[0]')
        print('this obsID was skipped', obsID)
        continue

    #for the titles on the plots to come up appropriately
    title_xmm = obsID
    overarching_obsIDs.append(obsID)

    print('########this obsID is going to be correlated###########', obsID)

    #trimming them
    time_xmm, xmm = zip(*[(t,y) for t, y in zip(time_xmm, xmm) if start <= t <= end])
    time_erosita, erosita = zip(*[(t,y) for t, y in zip(time_erosita, erosita) if start <= t <= end])

    #ensuring they all have the same length by cutting down to that length
    min_len = min(len(time_xmm), len(xmm),
              len(time_erosita), len(erosita))
    time_xmm = time_xmm[:min_len]
    xmm = xmm[:min_len]
    time_erosita = time_erosita[:min_len]
    erosita = erosita[:min_len]

    #as we will be using only one timestamp, need to ensure that the shift is then later accounted for again
    #(indeed there is only a shift as we ensured that the binning is otherwise the same)
    small_shift = time_xmm[0] - time_erosita[0]

    

    #first plot: just the two count rates
    if plot_wanted:
        fig_twin_1, ax1 = plt.subplots()

        ax1.plot(time_xmm, xmm, label="xmm", color = 'royalblue')
        ax1.set_ylabel('xmm count rate')
        ax1.set_xlabel('Time $s$')
        ax2 = ax1.twinx()
        ax2.plot(time_erosita, erosita, label="erosita", color = 'chocolate')
        ax2.set_ylabel('erosita count rate')
        fig_twin_1.suptitle(title_xmm)
        fig_twin_1.legend()

    #do not want to be doing correlation if there are less than x data points (here, chose 5)
    if len(time_erosita) <= 5:
        print('there were not enough datapoints to proceed with the correlation (<=5)')
        print('this obsID was skipped', obsID)
        if plot_wanted:
            plt.show()
        continue

    #########correlation for the whole dataset#########
    #print('the correlation for the whole dataset threshold was reached')
    delays_full_dataset, corr_full_dataset = cross_correlate(xmm, erosita, t_1=time_erosita, t_2=time_erosita)
    one_sig = np.std(corr_full_dataset)
    if plot_wanted:
        fig_main, axs_main = plt.subplots(1,2)
        for i in range(1,4):
            axs_main[0].hlines(i*one_sig, xmin=min(delays_full_dataset), xmax=max(delays_full_dataset), color="palevioletred")
            for i in range(1,4):
                i = -i
                axs_main[0].hlines(i*one_sig, xmin=min(delays_full_dataset), xmax=max(delays_full_dataset), color="palevioletred")

        axs_main[0].plot(delays_full_dataset, corr_full_dataset, color = 'black')
        axs_main[0].set_xlabel(r'Delays $\tau$')
        axs_main[0].set_ylabel(r'Correlation value $C(\tau)$')
        axs_main[0].set_title("Correlation plot for full data with multiple standard deviations as horizontal lines \n"+ title_xmm)

        #u.usual_plot(xl=r"Delays $\tau$", yl=r"Correlation value $C(\tau)$", title="Correlation plot for full data with multiple standard deviations as horizontal lines \n"+ title_xmm)
        axs_main[1].hist(corr_full_dataset, bins="auto", color = 'grey')
        axs_main[1].set_xlabel('Correlation')
        axs_main[1].set_ylabel('Frequency')
        axs_main[1].set_title("Distribution of $C(\tau)$ ")


    #u.usual_plot(xl=r"Correlation", yl="Frequency", title=r"Distribution of $C(\tau)$ ")
    #plt.show()

    #this provides the maximum delay for the whole set
    index_max, index_min = np.argmax(corr_full_dataset), np.argmin(corr_full_dataset)
    max_delay, min_delay = delays_full_dataset[index_max], delays_full_dataset[index_min]
    print('this is max', max_delay, 'this is min', min_delay)


    #want to make it such that it takes the arguments of all above 2sigma and hence the delays
    
    high_correlation_idx = [np.where(corr_full_dataset == c) for c in corr_full_dataset if c >= 2*one_sig]
    high_correlation_delays = [delays_full_dataset[i] for i in high_correlation_idx]
    print('these are the high correlation delays', high_correlation_delays)
    print('this is the small shift between the xmm and erosita time series (time_xmm[0] - time_erosita[0])', small_shift)
    time_erosita = np.array(time_erosita)
    optimal_delay = max_delay + small_shift #such that we can plot on erosita time - if plotting both on xmm timestamps then just optimal_delay = max +delay
    #sub_optimal_delay = second_max + small_shift

    
    if corr_full_dataset[index_max] > 2*one_sig:
        more_than_2stds.append(optimal_delay)

    #if corr_full_dataset[index_max] >3*one_sig:
    #    more_than_3stds.append([optimal_delay, obsID])

    #instead trying to do it by including every delay above the 3 sigma, since each may be a different flare: allowing for more than one per flare
    for delay,index in zip(high_correlation_delays,high_correlation_idx):                                           
        delay = int(delay[0]) #because for some reason the delays came as arrays of length 0
        if corr_full_dataset[index] > 3*one_sig:
            more_than_3stds.append([delay,obsID])



    #delays from the overarching correlation
    overarching_delays = [h + small_shift for h in high_correlation_delays]

    if plot_wanted and len(overarching_delays) != 0 and len(overarching_delays) != 1:
        fig, axs = plt.subplots(nrows = len(overarching_delays), ncols = 1)
        for i in range(len(overarching_delays)):
            print('this is i', i)
            axs[i].plot(time_erosita, xmm, label = 'xmm', color = 'royalblue')
            axs_2 = axs[i].twinx()
            axs_2.plot(time_erosita - overarching_delays[i], erosita, label= 'erosita', color = 'chocolate')
            fig.suptitle("All delays with > 2sigma correlation value, \n from correlating the dataset as a whole")
            fig.legend()
            print('plot with this delay should have been pdt', overarching_delays[i])

    #hashed it since removed the second max specifically to have the automatic array with all those above 2 sigma instead
    if plot_wanted:
        fig_shifted_opt, axs_shifted_opt_1 = plt.subplots()
        axs_shifted_opt_1.plot(time_erosita, xmm, label = 'xmm', color = 'royalblue')
        axs_shifted_opt_1.set_ylabel('xmm count rate')
        axs_shifted_opt_2 = axs_shifted_opt_1.twinx()
        axs_shifted_opt_2.plot(time_erosita - optimal_delay,erosita, label = 'erosita', color = 'chocolate')
        axs_shifted_opt_2.set_ylabel('erosita count rate')
        fig_shifted_opt.suptitle("SHIFTED time series \n" + title_xmm+ "\n with max shift from full data set")
        fig_shifted_opt.legend()

    print('#####################################################################\n this is optimal_delay in hours, mins, secs', int(optimal_delay//3600), 'h', int((optimal_delay%3600)//60),'min', optimal_delay %60, 's\n###################################################################################')


    
    overarching_optimal_delays.append(optimal_delay)

    #print('this is sub_optimal_delay in hours, mins, secs', int(sub_optimal_delay/3600), 'h', int((sub_optimal_delay%3600)//60),'min', sub_optimal_delay %60, 's')






   





                                                                                                                            

    if plot_wanted:
        plt.show()
print('############################################################################################################################################')
print('these are all the overarching delays (optimised)',overarching_optimal_delays)
print('these are all delays with corr > 2 sigma', more_than_2stds)
print('these are all the delays with corr > 3 sigma', more_than_3stds)
print('these are the corresponding ObsIDs', overarching_obsIDs)

with open("correlation_ultimate_output_FINAL_with_degeneracy.json", "w") as f:
        json.dump(overarching_obsIDs, f, indent = 2)
        json.dump(overarching_optimal_delays,f,indent = 2)
        json.dump(more_than_2stds, f, indent = 2)
        json.dump(more_than_3stds,f,indent = 2)
        print('the monster list was dumped: contains in order overarching_obsIDs, optimal_delays, more_than_2stds delays')
