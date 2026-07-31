import numpy as np
import json
import matplotlib.pyplot as plt

#initialising
large_list = []
plot_wanted = False
plot_poster = False

#getting the data from all_eek.json
with open("all_eek.json", "r") as f:
    large_list = json.load(f)
    large_list_mod = [large_list[0],0,0,0]

    #so this large list has, in order
    #obs IDs, all delays, delays with corr > 2std, delays with corr > 3std

    #print('do my elements have the same length??', len(large_list[0]), len(large_list[1]), len(large_list[2]), len(large_list[3]))

    #filtering out irrelevant numbers
    large_list_mod[1] = [n for n in large_list[1] if n < 10000]
    large_list_mod[1] = [n for n in large_list[1]]
    large_list_mod[2] = [n for n in large_list[2] if n < 10000]
    large_list_mod[3] = [n for n in large_list[3] if n < 10000]
    small_delays_with_3std = [n for n in large_list[3] if n < 120]
    print('these are the len of the really small delays with 3std', len(small_delays_with_3std))
    print('these are the lengths of all', len(large_list_mod[1]), 'and above 3 sigma', len(large_list_mod[3]))
    print('this is the ratio', len(large_list_mod[3])/len(large_list_mod[1]))


if plot_poster:
    bin_nb = 14
    fig, ax = plt.subplots(figsize=(4.5, 6))
    ax.hist(large_list_mod[3], bins=bin_nb, color='#649B90', edgecolor='white', linewidth=0.5)

    #ax.set_xlabel(r'Optimal delay $\tau$ with $C(\tau) > 3\sigma$', fontsize=12)
    ax.set_xlabel(r'$\tau$', fontsize=12)
    ax.set_ylabel('Number of correlated events', fontsize=12)

    ax.tick_params(direction='in', length=6, width=1.5, labelsize=12, which='both')
    ax.tick_params(length=3, which='minor')
    ax.minorticks_on()
    ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=6, prune='both'))

    for spine in ax.spines.values():
        spine.set_linewidth(1.2)

    plt.tight_layout()
    plt.savefig('3sigma_histogram.png', dpi=150)
    plt.show()




#plotting histograms: binning how many obsIDs had a certain correlation length
if plot_wanted:
    bin_nb = 12
    fig, axs = plt.subplots(nrows = 1, ncols = 4)
    axs[0].hist(large_list_mod[1], bins = bin_nb, color = 'orchid')
    axs[0].set_xlabel('Optimal delay for each obsID')
    axs[1].hist(large_list_mod[2], bins = bin_nb, color = 'orchid')
    axs[1].set_xlabel(r'Optimal delays with corr > $2\sigma$')
    axs[2].hist(large_list_mod[3], bins = bin_nb, color = 'orchid')
    axs[2].set_xlabel(r'Optimal delays with corr > $3\sigma$')
    axs[3].hist(small_delays_with_3std, color = 'darkorchid')
    axs[3].set_xlabel(r'Optimal small delays with corr > $3\sigma$')
    plt.show()


#next step: go and investigate the super long ones and how that offset even works. Maybe try setting a max to the offset if it looks irrelevant.
#go and look at the correlation plots of some of the >3std


#creating a list and adding the eROday before and after
#for it to work: just need to add a useless value at the end as it will not analyse it

#all the eROdays with enhancements
list_2 = [43205, 43211, 43216, 43222, 43230, 43231, 43233, 43234, 43253, 43259, 43265, 43267, 43268, 43269, 43272, 43277, 43283, 43289, 43291, 43294, 43296, 43297, 43298, 43300, 43302, 43303, 43318, 43319, 43322, 43334, 43338, 43341, 43360, 43386, 43399, 43403, 43415, 43421, 43424, 43425, 43426, 43427, 43428, 43431, 43432, 43433, 43434, 43435, 43436, 43437, 43440, 43443, 43444, 43445, 43446, 43447, 43448, 43449, 43450, 43454, 43459, 43463, 43464, 43467, 43546, 43578, 43582, 43588, 43590, 43591, 43593, 43595, 43597, 43598, 43631, 43645, 43684, 43703, 43705, 43710, 43733, 43835, 43853, 43869, 43916, 44007, 44009, 44011, 44012, 44052, 44054, 44057, 44093, 44124, 44144, 44243, 44249, 44253, 44254, 44290, 44291, 44302, 44307, 44312, 44338, 44343, 44345, 44371, 44372, 44378, 44383, 44386, 44402, 44410, 44423, 44424, 44425, 44439, 44442, 44446, 44447, 44448, 44449, 44450, 44452, 44453, 44454, 44458, 44462, 44469, 44511, 44520, 44521, 44533, 44561, 44567, 44627, 44673, 44674, 44727, 44780, 44956, 44957, 44963, 44987, 44988, 44991, 44998, 45012, 45122, 45124, 45126, 45133, 45137, 45139, 45144, 45227, 45270, 45275, 45277, 45279, 45292, 45293, 45294, 45296, 45297, 45298, 45397, 45415, 45433, 45434, 45436, 45439, 45446, 45449, 45450, 45457, 45460, 45463, 45464, 45470, 45471, 45475, 45479, 45492, 45515, 45516, 45521, 45524, 45525, 45571, 45582, 45613, 45615, 45618, 45619, 45622, 45629, 45630, 45641, 45646, 45647, 45659, 45664, 45683, 45684, 45695, 45701, 45722, 45750, 45756, 45763, 45779, 45810, 45827, 45832, 45833, 45840, 45845, 45851, 45855, 45858, 45861, 45867, 45869, 45873, 45879, 45885, 45888, 45890, 45891, 45892, 45893, 45894, 45897, 45898, 45899, 45901, 45902, 45903, 45904, 45941, 45948, 45976, 45983, 46001, 46058, 46060, 46076, 46077, 46096, 46160, 46161, 46163, 46190, 46231, 46251, 46252, 46253, 46254, 46258, 46259, 46288, 46291, 46313, 46314, 46326, 46329, 46332, 46333, 46334, 46335, 46336, 46337, 46341, 46342, 46343, 46346, 46347, 46349, 46350, 46351, 46358, 46372, 46415, 46417, 46418, 46420, 46422, 46423, 46427, 46451, 46452, 46459, 46460, 46462, 46464, 46465, 46480, 46492, 46493, 46497, 46498, 46500, 46501, 46502, 46503, 46504, 46505, 46507, 46509, 46511, 46512, 46516, 46518, 46520, 46541, 46558, 46576, 46578, 46590, 46607, 46610, 46633, 46635, 46637, 46665, 46666, 46667, 46669, 46671, 46672, 46673, 46674, 46675, 46677, 46678, 46679, 46680, 46681, 46682, 46684, 46690, 46697, 46698, 46699, 46700, 46701, 46705, 46711, 46712, 46713, 46717, 46727, 46751, 46788, 46807, 46808, 46809, 46813, 46814, 46816, 46819, 46820, 46835, 46848, 46857, 46861, 46865, 46869, 46870, 46871, 46879, 46881, 46882, 46890, 46892, 46896, 46899, 46900, 46901, 46902, 46903, 46917, 46921, 46923, 46924, 46925, 46926, 46929, 46934, 46935, 46946, 46947, 46950, 46960, 46971, 46975, 46977, 47024, 47028, 47030, 47110, 47111, 47118, 47127, 47133, 47166, 47167, 47168, 47169, 47176, 47177, 47178, 47219, 47241, 47242, 47279, 47281, 47289, 47316, 47340, 47343, 47358, 47359, 47382, 47389, 47409, 47413, 47422, 47423, 47424, 47425, 47426, 47427, 47429, 47431, 47433, 47436, 47458, 47463, 47476, 47490, 47506, 47518, 47525, 47526, 47537, 47539, 47543, 47565, 47567, 47568, 47572, 47577, 47578, 47579, 47580, 47582, 47586, 47587, 47612, 47613, 47616, 47618, 47620, 47621, 47634, 47639, 47642, 47644, 47651, 47653, 47656, 47658, 47659, 47660, 47661, 47663, 47664, 47665, 47666, 47667, 47668, 47670, 47680, 47693, 47694, 47715, 47724, 47725, 47726, 47727, 47730, 47731, 47732, 47733, 47734, 47736, 47740, 47741, 47746, 47751, 47777, 47783, 47803, 47808, 47816, 47819, 47820, 47829, 47830, 47831, 47834, 47838, 47840, 47841, 47842, 47843, 47845, 47846, 47847, 47849, 47851, 47852, 47853, 47854, 47855, 47856, 47857, 47858, 47859, 47860, 47861, 47862, 47865, 47867, 47872, 47873, 47874, 47878, 47879, 47880, 47882, 47883, 47884, 47887, 47903, 47905, 47913, 47924, 47938, 47941, 47942, 47943, 47944, 47945, 47946, 47947, 47948, 47950, 47952, 47956, 47957, 47958, 47959, 47968, 47969, 47971, 47972, 47974, 47975, 47983, 47984, 47999, 48001, 48009, 48011, 48015, 48017, 48020, 48021, 48033, 48035, 48037, 48043, 48057, 48058, 48076, 48087, 48088, 48119, 48127, 48144, 48146, 48148, 48151, 48152, 48153, 48154, 48178, 48197, 48206, 48207, 48208, 48209, 48213, 48227, 48234, 48261, 48264, 48267, 48268, 48269, 48272, 48296, 48300, 48303, 48304, 48305, 48308, 48309, 48310, 48311, 48312, 48317, 48322, 48324, 48328, 48330, 48331, 48332, 48333, 48334, 48335, 48336, 48337, 48339, 48345, 48351, 48356, 48357, 48358, 48359, 48364, 48370, 48378, 48389, 48391, 48398, 48399, 48400, 48402, 48407, 48408, 48409, 48410, 48411, 48414, 48415, 48417, 48418, 48423, 48425, 48437, 48458, 48462, 48469, 48481, 48486, 48498, 48509, 48510, 48511, 48512, 48513, 48514, 48518, 48523, 48526, 48539, 0]
'''
def create_large_list():
    large_list = []
    for i in range(len(list_2) -1):
        before = 0
        of = list_2[i]
        after = 0
        if list_2[i-1] != list_2[i]-1:
            large_list.append(list_2[i]-1)
        large_list.append(list_2[i])
        if list_2[i+1] != list_2[i]+1:
            large_list.append(list_2[i]+1)
    print(list(set(large_list)))
'''



#have a file with all the obsIDs coupled to their durations.
#Need to match all the all.json output to their duration
with open("obsID_durations.json", "r") as f:
    durations = json.load(f)
    print('this is len durations', len(durations))
durations_of_enhancement_obsIDs = []
for element in large_list[0]:
    index = [i for i, obs in enumerate(durations) if element == obs[0]]
    corresponding_element =durations[index[0]]
    corresponding_duration = corresponding_element[1]
    durations_of_enhancement_obsIDs.append(corresponding_duration)

#print('do my elements have the same length??', len(large_list[0]), len(large_list[1]), len(large_list[2]), len(large_list[3]))
#print(np.max(durations_of_enhancement_obsIDs))
#print('checking the correpsonding durations are indeed the same lenth as the actual enhanced and correlated obsIDs')
#print(len(large_list[0]),len(large_list[1]),len(durations_of_enhancement_obsIDs))

'''
if plot_wanted:
    fig, axs = plt.subplots()
    axs.scatter(durations_of_enhancement_obsIDs,large_list[1])
    plt.xlabel('Duration of the obsID (s)')
    plt.ylabel('Optimal delays (s)')
    plt.ylim(ymin = 0, ymax = 10000)
    plt.show()
    '''


#coupling the obsID name to whether or not eRO and XMM were on the same side of the ecliptic
with open("positions_eek.json", "r") as f:
    positions = json.load(f)
    #print('this is the length of positions', len(positions[0]), len(positions[1]))
    coupled_truths = []
    for truth, element in zip(positions[0], positions[1]):
        single_coupled_truth = [truth, element[0]]
        coupled_truths.append(single_coupled_truth)
    print('this is the length of coupled truth', len(coupled_truths))
    print('this is the first element of coupled truth', coupled_truths[0])
    #print(coupled_truths)

def match_reduced_to_obsID(total_coupled, reduced_huh): #i.e. the full obsID list, the full delays list and the reduced delay list
    print(total_coupled[1])
    obsIDs_huh = [truth[0] for truth in total_coupled]
    print(obsIDs_huh[0])
    optimals_huh = [truth[1] for truth in total_coupled]
    #we can analogously do with the full truths list, the full obsIDs list and the reduced obsID list
    queen = []
    for freddie_mercury in reduced_huh: #because it has >2std correlation
        index = [i for i, delay in enumerate(optimals_huh) if freddie_mercury == delay]
        corresponding_obsID = obsIDs_huh[index[0]]
        queen.append(corresponding_obsID)
    return queen #here queen will be the reduced truths

reduced_truths = match_reduced_to_obsID(coupled_truths, large_list[0])

#plotting the obsID length vs the correltion length
# andcolour coded for above or below the ecliptic (i.e. same side or not)
if plot_wanted:
    fig, axs = plt.subplots()
    for truth, delay,duration in zip(reduced_truths, large_list[1],durations_of_enhancement_obsIDs):
        if truth == True:
            plt.scatter(duration,delay, color = 'green')
        elif truth == False:
            plt.scatter(duration, delay, color = 'red')
        else:
            continue

    plt.xlabel('Duration of the obsID (s)')
    plt.ylabel('Optimal delays (s)')
    plt.ylim(ymin = 0, ymax = 10000)
    plt.show()


#new part added from Claude
plot_wanted_SPIE = True

#last iteration of the code which yields every observationID with C(tau) > 3 sigma with its correspoding optimal delay

#only keep the last one above 3 sigma
new_list_1 = [[469.7438780523481, '0831801401'], [147.06417721285482, '0842110101'], [321.9102621489558, '0844100701'], [361.6178332169851, '0810850501'], [37.53085434436798, '0842591001'], [3304.667142083079, '0844930101'], [1795.6563273056418, '0840841001'], [670.8612038950835, '0841180401'], [19660.1411125809, '0844100101'], [3953.095799675639, '0842570101'], [17054.877530279824, '0844020101'], [20755.0865184573, '0851182001'], [8191.507839490668, '9364100002'], [1644.7681964244998, '9364100003'], [12630.611022675008, '0851181401'], [346.80223563171563, '0823810301'], [7444.038026434027, '0843020801'], [10.130776047706604, '0840131201'], [7096.184145745525, '0843830401'], [96.34067225456238, '0840580401'], [24766.999936695403, '0827230701'], [102.10020077228546, '0827230801'], [1488.5667508052243, '9370900013'], [252.74813189958547, '0840133901'], [24113.63716459842, '0840750201'], [2736.595282234397, '0840842901'], [15588.699475517407, '0840211101'], [7562.78527709282, '0852060301'], [4403.447498827702, '0843440101'], [937.6761287263795, '0841570401'], [1313.6979128360747, '0841480201'], [519.4839675520597, '0854591201'], [5489.592575939789, '0842730101'], [1032.6853335600981, '0827251201'], [2471.824180462914, '0863560601'], [887.6451124668121, '0863420401'], [71.12647271156311, '0827231801'], [1852.1643586024159, '0860370201'], [285.2964832481188, '0861910301'], [1495.615673794303, '0865040401'], [10596.253605353206, '0811023801'], [56349.97743117468, '0862220101'], [182.76522057464248, '0827050601'], [16585.501842251033, '0871010101'], [352.62488877432685, '0860140101'], [110.32872557640076, '0827060701'], [5130.237789260016, '0862470801'], [1876.2986605809283, '0862470901'], [7521.814323504766, '0862471001'], [2081.4794623609228, '0862471101'], [885.4685489384402, '0862220301'], [12455.921470993755, '0864622301'], [26872.14554106848, '0864621401'], [4156.620046193517, '0870790201'], [12507.503296969897, '0861080201'], [1437.3542125853044, '0823594201'], [95.38752448558807, '0863560501'], [5813.508774938869, '0860303601'], [148.92954139344326, '0864621501'], [1564.5112477949688, '0870930301'], [33.723228216171265, '0871590201'], [208.34249366469038, '0824320701'], [1672.0194052815436, '9383000004'], [835.7631548656567, '0864090301'], [1518.7760949307476, '0862641601'], [12428.146777707701, '0864622001'], [9348.606408913158, '0861680201'], [1432.1057397270204, '0865011001'], [6188.30733389968, '9384200002'], [8182.423879463117, '0864410101'], [3602.158392577271, '0863400601'], [2649.9663751484986, '0827251501'], [40116.810283535204, '0827251701'], [569.4891147707024, '9384700003']]

#keep every delay above 3 sigma: allow for degeneracy
new_list_2 = [[367, '0831801401'],  [120, '0842110101'], [362, '0842110101'], [482, '0842110101'],  [120, '0844100701'], [241, '0844100701'], [362, '0844100701'], [482, '0844100701'], [242, '0810850501'], [121, '0842591001'], [242, '0842591001'], [363, '0842591001'], [3042, '0844930101'], [3163, '0844930101'], [3285, '0844930101'], [3407, '0844930101'], [1696, '0840841001'], [485, '0841180401'], [606, '0841180401'], [18709, '0844100101'], [18829, '0844100101'], [19554, '0844100101'], [3870, '0842570101'], [14635, '0842570101'], [16986, '0844020101'], [22648, '0844020101'], [20591, '0851182001'], [20712, '0851182001'], [20832, '0851182001'], [36366, '0851182001'], [36486, '0851182001'], [36607, '0851182001'], [36727, '0851182001'], [36848, '0851182001'], [362, '9364100002'], [725, '9364100002'], [8096, '9364100002'], [1588, '9364100003'], [1711, '9364100003'], [12173, '0851181401'], [12293, '0851181401'], [12534, '0851181401'], [12655, '0851181401'], [12896, '0851181401'], [246, '0823810301'], [1807, '0843020801'], [1928, '0843020801'], [2048, '0843020801'], [2169, '0843020801'], [2289, '0843020801'], [2410, '0843020801'], [2530, '0843020801'], [6989, '0843020801'], [7109, '0843020801'], [7230, '0843020801'], [7350, '0843020801'], [7471, '0843020801'], [7591, '0843020801'], [7712, '0843020801'], [7832, '0843020801'],  [6154, '0843830401'], [6275, '0843830401'], [6758, '0843830401'], [6878, '0843830401'], [6999, '0843830401'], [7120, '0843830401'], [7240, '0843830401'],  [120, '0840580401'], [241, '0840580401'], [7825, '0827230701'], [7945, '0827230701'], [24679, '0827230701'],  [120, '0827230801'], [240, '0827230801'], [361, '0827230801'], [241, '9370900013'], [362, '9370900013'], [845, '9370900013'], [1448, '9370900013'], [1569, '9370900013'], [1810, '9370900013'], [241, '0840133901'], [6749, '0840750201'], [7834, '0840750201'], [24106, '0840750201'], [2697, '0840842901'], [15559, '0840211101'], [6279, '0852060301'], [6400, '0852060301'], [7366, '0852060301'], [7487, '0852060301'], [7608, '0852060301'], [7728, '0852060301'], [4338, '0843440101'], [25788, '0843440101'], [607, '0841570401'], [849, '0841570401'], [966, '0841480201'], [1208, '0841480201'], [120, '0854591201'], [482, '0854591201'], [965, '0854591201'], [3377, '0854591201'], [13993, '0854591201'], [5426, '0842730101'], [12058, '0842730101'], [728, '0827251201'], [849, '0827251201'], [970, '0827251201'], [1092, '0827251201'], [1213, '0827251201'], [2295, '0863560601'], [2416, '0863560601'], [3624, '0863560601'], [3745, '0863560601'], [853, '0863420401'], [16278, '0827231801'], [1813, '0860370201'], [1934, '0860370201'], [2054, '0860370201'], [2296, '0860370201'], [2417, '0860370201'], [3626, '0860370201'], [241, '0861910301'], [6398, '0861910301'], [1459, '0865040401'], [10522, '0811023801'], [10643, '0811023801'], [56252, '0862220101'], [121, '0827050601'], [605, '0827050601'], [16167, '0871010101'], [16288, '0871010101'], [16529, '0871010101'],  [120, '0860140101'], [241, '0860140101'], [362, '0860140101'],  [604, '0862470801'], [4955, '0862470801'], [5076, '0862470801'], [5197, '0862470801'], [5922, '0862470801'], [1849, '0862470901'], [7488, '0862471001'], [7609, '0862471001'], [2052, '0862471101'], [845, '0862220301'], [12365, '0864622301'], [10966, '0864621401'], [11086, '0864621401'], [25306, '0864621401'], [25426, '0864621401'], [25547, '0864621401'], [26752, '0864621401'], [26872, '0864621401'],  [3135, '0870790201'], [3255, '0870790201'], [3376, '0870790201'], [3496, '0870790201'], [3617, '0870790201'], [3738, '0870790201'], [3858, '0870790201'], [3979, '0870790201'], [4099, '0870790201'], [4220, '0870790201'], [4340, '0870790201'], [4461, '0870790201'], [4582, '0870790201'], [4702, '0870790201'], [4823, '0870790201'], [12434, '0861080201'], [1085, '0823594201'], [1206, '0823594201'], [1326, '0823594201'], [1567, '0823594201'], [1688, '0823594201'], [1809, '0823594201'], [1929, '0823594201'], [2050, '0823594201'], [120, '0863560501'], [241, '0863560501'], [361, '0863560501'], [482, '0863560501'], [723, '0863560501'], [5793, '0860303601'], [121, '0864621501'], [1456, '0870930301'], [2306, '0870930301'],  [131, '0871590201'],  [120, '0824320701'], [241, '0824320701'], [361, '0824320701'], [482, '0824320701'], [602, '0824320701'], [723, '0824320701'], [843, '0824320701'], [964, '0824320701'], [1084, '0824320701'], [1205, '0824320701'], [1325, '0824320701'], [1446, '0824320701'], [1575, '9383000004'], [602, '0864090301'], [723, '0864090301'], [1084, '0864090301'], [2289, '0864090301'], [2530, '0864090301'], [2651, '0864090301'], [2892, '0864090301'], [3253, '0864090301'], [3615, '0864090301'], [3736, '0864090301'], [3856, '0864090301'], [847, '0862641601'], [1452, '0862641601'], [12426, '0864622001'], [9061, '0861680201'], [9182, '0861680201'], [9303, '0861680201'], [9424, '0861680201'], [9544, '0861680201'], [1348, '0865011001'], [361, '9384200002'], [5901, '9384200002'], [6021, '9384200002'], [6141, '9384200002'], [6262, '9384200002'], [6382, '9384200002'], [8095, '0864410101'], [3505, '0863400601'], [2533, '0827251501'], [39636, '0827251701'], [39757, '0827251701'], [39877, '0827251701'], [39998, '0827251701'], [40119, '0827251701'], [40240, '0827251701'], [40361, '0827251701'], [40482, '0827251701'], [40602, '0827251701'], [41327, '0827251701'], [482, '9384700003'], [2049, '9384700003'], [2170, '9384700003'], [2411, '9384700003'], [3135, '9384700003']]


#delays above 3 sigma (I forget whether degeneracy or not, will check) - INCLUDING BACKWARDS
#to hopefully do
positions_truth_FINAL =  [[True, '0831801401'], [True, '0810841401'], [True, '0831801501'], [True, '0851181301'], [True, '0784303201'], [True, '0842110101'], [True, '0843151101'], [True, '0843150501'], [True, '9362600002'], [True, '0844100701'], [True, '0840132201'], [True, '0852180501'], [True, '0850780201'], [True, '0841180101'], [True, '0840841101'], [True, '0841320201'], [True, '0810850501'], [True, '0853980401'], [True, '0853980301'], [True, '0840720101'], [True, '0840210301'], [True, '0842362401'], [True, '0842362101'], [True, '0842591001'], [True, '9363000003'], [True, '0844930101'], [True, '0840841001'], [True, '9363100004'], [True, '0841510101'], [True, '0841180401'], [True, '0853220301'], [True, '0853000201'], [True, '0844100101'], [True, '0841920101'], [True, '0840910801'], [True, '0842570101'], [True, '0844020101'], [True, '9364000003'], [True, '0852980501'], [True, '0851182001'], [True, '0851181601'], [True, '9364100002'], [True, '9364100003'], [True, '0853980201'], [True, '0851181401'], [True, '0823810301'], [True, '0843151001'], [True, '0853980601'], [True, '0843020801'], [True, '0841230401'], [True, '9364300005'], [True, '0804070101'], [True, '0810890201'], [True, '9365300003'], [True, '0853981101'], [True, '0853980901'], [True, '0853981001'], [True, '0840131201'], [True, '0843830401'], [True, '0844860201'], [False, '0853782501'], [False, '0844860501'], [False, '0844860601'], [False, '0840580401'], [False, '0853230301'], [False, '0827221901'], [False, '0827230701'], [False, '0827230801'], [False, '0840740101'], [False, '0870590101'], [False, '9370900013'], [True, '0853790101'], [True, '0841450101'], [True, '0842550101'], [True, '0840490101'], [True, '0840133901'], [True, '9371700002'], [True, '0865600201'], [True, '0844970601'], [True, '9371900002'], [True, '0840750201'], [True, '0810841501'], [True, '0840842901'], [True, '0840211101'], [True, '0841951701'], [True, '0840843001'], [True, '0852060301'], [True, '0841660301'], [True, '0840841201'], [True, '0843440101'], [True, '0840841301'], [True, '9372600004'], [True, '0841570401'], [True, '0841480201'], [True, '0810871301'], [True, '0840550101'], [True, '0854591201'], [True, '0842730101'], [True, '0810863001'], [True, '0827251201'], [True, '0863560601'], [True, '0864330201'], [True, '0863810301'], [True, '0810821501'], [True, '0863420401'], [True, '0864560101'], [True, '0862731101'], [True, '0870810101'], [False, '0827231801'], [False, '0860370201'], [False, '0861910301'], [False, '0872390101'], [False, '0811023901'], [False, '0865040401'], [False, '0811023801'], [False, '0862220101'], [False, '9379600003'], [False, '0827050601'], [False, '0862470401'], [False, '0860303001'], [False, '0871010101'], [False, '0861111201'], [False, '0860140101'], [False, '0827060701'], [False, '0862470801'], [False, '0862470901'], [False, '0862471001'], [True, '0862471101'], [True, '0862220301'], [True, '0810850601'], [True, '0861310101'], [True, '0864622301'], [True, '0864621401'], [True, '0871191301'], [True, '0870790201'], [True, '0861080201'], [True, '9381500004'], [True, '0823594201'], [True, '0863560501'], [True, '0860303601'], [True, '0872990401'], [True, '0864621501'], [True, '0870870601'], [True, '0870870201'], [True, '0871191401'], [True, '0862640201'], [True, '0870930701'], [True, '0870930301'], [True, '9382600002'], [True, '0861950201'], [True, '0871590201'], [True, '0870930401'], [True, '0824320701'], [True, '9383000004'], [True, '0865450401'], [True, '0864090301'], [True, '9383500002'], [True, '0864622001'], [True, '0865010101'], [True, '0862920201'], [True, '0861680201'], [True, '0865011001'], [True, '9384200002'], [True, '0862840101'], [True, '0865011201'], [True, '9384300002'], [True, '0827350401'], [True, '0865011101'], [True, '0864410101'], [True, '0863400601'], [True, '0865400201'], [True, '0827251501'], [True, '0827360501'], [True, '0827251601'], [True, '9384600003'], [True, '0872391201'], [True, '0861582401'], [True, '0871020101'], [False, '9384700002'], [False, '0827251701'], [False, '9384700003'], [False, '0870190501'], [False, '0827321301'], [False, '0860800301'], [False, '0863780201'], [False, '0862730101'], [False, '0865350201'], [False, '0827361101'], [False, '0862150101'], [False, '0864052301'], [False, '0863230301'], [False, '0872391401'], [False, '0860530201'], [False, '0863960101'], [False, '0871390501'], [False, '0864081201'], [False, '0840310201'], [False, '0862950801'], [False, '0884220201'], [False, '9388300002'], [False, '0862900501'], [False, '0861582001'], [False, '0811012701'], [False, '0870850201'], [False, '0811024001'], [False, '0862980101'], [False, '0852600201'], [False, '0862900201'], [False, '0862950401'], [False, '0865140801'], [False, '0862900301'], [False, '0862770701'], [False, '0860302601'], [False, '0852600301'], [False, '0886010501'], [False, '0864621301'], [False, '0861171301'], [False, '0864550501'], [False, '0864530101'], [False, '0861171501'], [False, '0861171201'], [False, '0853210501'], [False, '0886010201'], [False, '0861172001'], [False, '0861171901'], [False, '0864550601'], [False, '0860620401'], [False, '0861171601'], [False, '9390300004'], [True, '0861260101'], [True, '0861171101'], [True, '0864340101'], [True, '0864440201'], [True, '0870940401'], [True, '0810871401'], [True, '0872392801'], [True, '0882340601'], [True, '0863400301'], [True, '0882870201'], [True, '0870870801'], [True, '0864050601'], [True, '0864052501'], [True, '0865011301'], [True, '9391500002'], [True, '9391500003'], [True, '0863090101'], [True, '0884710101'], [True, '0865380201'], [True, '0882110401'], [True, '0865011601'], [True, '0872392901'], [True, '0865011701'], [True, '0883950101'], [True, '0810811601'], [True, '0882160401'], [True, '0882720301'], [True, '0864430101'], [True, '0864430201'], [True, '0870830401'], [True, '0881990101'], [True, '0881561101'], [True, '0881560601'], [True, '0871590701'], [True, '0883780101'], [True, '9393200007'], [True, '0880030501'], [True, '0865012801'], [True, '0865011801'], [True, '9393300002'], [True, '0884960101'], [True, '0882110501'], [True, '0864010101'], [True, '0890400301'], [True, '0881071701'], [True, '0882110601'], [True, '0884370101'], [True, '9394100002'], [True, '0842760201'], [True, '0880710301'], [True, '9394800002'], [True, '0880000301'], [True, '9394800003'], [True, '0884970101'], [False, '0880000401'], [False, '0885100101'], [False, '0884190401'], [False, '0891020101'], [False, '0882340701'], [False, '9395900003'], [False, '0882260101'], [False, '0880580601'], [False, '0884991701'], [False, '0883800101'], [False, '0890660101'], [False, '0862091101'], [False, '0870840101'], [False, '0862090501'], [False, '9397100005'], [False, '0882060701'], [False, '0862091001'], [False, '0882870601'], [False, '0891801101'], [False, '0862770501'], [False, '0882061401'], [False, '0882060901'], [False, '0882870501'], [False, '9397500003'], [False, '0882640301'], [False, '0880760101'], [False, '0891801501'], [False, '0861470401'], [False, '0861470201'], [False, '0882650101'], [False, '0886020401'], [False, '0880280801'], [False, '0882061201'], [False, '0883500201'], [False, '0883500101'], [False, '0890420101'], [False, '0884860101'], [False, '0886041101'], [False, '0883460301'], [False, '9398600003'], [False, '0881680301'], [False, '0841930601'], [False, '0841930501'], [False, '0881840401'], [False, '0880280401'], [False, '0880280301'], [False, '0891801801'], [False, '0891801701'], [False, '0884080101'], [False, '0884050101'], [False, '0882031101'], [False, '0882030501'], [False, '0890450101'], [False, '0880281201'], [False, '0880280101'], [False, '0885090201'], [False, '0885090101'], [False, '0890440101'], [False, '0880860601'], [False, '0880031301'], [False, '0810841901'], [False, '0881680901'], [False, '0881681101'], [False, '0880810301'], [True, '0880810101'], [True, '0880810401'], [True, '0891010101'], [True, '0891010201'], [True, '0882340801'], [True, '0810871501'], [True, '9400900002'], [True, '0810880801'], [True, '0810880901'], [True, '0882721001'], [True, '0882720201'], [True, '9401000002'], [True, '0890810201'], [True, '9401000003'], [True, '0884960201'], [True, '0891802601'], [True, '9401000005'], [True, '0891802701'], [True, '0871592201'], [True, '0871591301'], [True, '9401100006'], [True, '0881680201'], [True, '0810243201'], [True, '0891070201'], [True, '0891070101'], [True, '9401200002'], [True, '0882030701'], [True, '0885240201'], [True, '0885240101'], [True, '0810890501'], [True, '0892420201'], [True, '9401500003'], [True, '0881350101'], [True, '0881210101'], [True, '0883090601'], [True, '0883090301'], [True, '9401800002'], [True, '0891802101'], [True, '0881210201'], [True, '0891801301'], [True, '0883780201'], [True, '0862860301'], [True, '0883210401'], [True, '0881070301'], [True, '0883550101'], [True, '0883550201'], [True, '0891803001'], [True, '0893400101'], [True, '0882721101'], [True, '0882720501'], [True, '9402400002'], [True, '0872610401'], [True, '0872610701'], [True, '9402600002'], [True, '0810850701'], [True, '0891803801'], [True, '0883040201'], [True, '0881900301'], [True, '0884470201'], [True, '0891804201'], [True, '0881072801'], [True, '0893810101'], [True, '0891802401'], [True, '0884993801'], [True, '9403900002'], [True, '0872610901'], [True, '0872610601'], [True, '0893810201'], [True, '0881900801'], [True, '0884992901'], [True, '0883550301'], [True, '0890640201'], [True, '0880000701'], [True, '0881071301'], [True, '0884120101'], [True, '0891400201'], [True, '0884120201'], [True, '0884740201'], [True, '0881190101'], [True, '0893810401'], [True, '0885260101'], [False, '0890700201'], [False, '0882050401'], [False, '9405100003'], [False, '0881190301'], [False, '9405100004'], [False, '0881990201'], [False, '0880860201'], [False, '0880001201'], [False, '9405200004'], [False, '0870840201'], [False, '0870860401'], [False, '0882590901'], [False, '9405300002'], [False, '0884630201'], [False, '0882650201'], [False, '0882650301'], [False, '0840140401'], [False, '9405600002'], [False, '0893810601'], [False, '0892420401'], [False, '0840140701'], [False, '0840140901'], [False, '9405700002'], [False, '0880900101'], [False, '0893810801'], [False, '0840140801'], [False, '0840140601'], [False, '0890700301'], [False, '0885260201'], [False, '0882130301'], [False, '9406300003'], [False, '0881770201'], [False, '0882650501'], [False, '0884120301'], [False, '9406500003'], [False, '0890620101'], [False, '0880582201'], [False, '0810401301'], [False, '0885110101'], [False, '0881630601']]



with open("obsID_durations.json", "r") as f:
    durations = json.load(f)

duration_dict = {obs[0]: obs[1] for obs in durations}
truth_dict = {entry[1]: entry[0] for entry in positions_truth_FINAL}

colour_map = {True: 'green', False: 'red'}

red, green, grey = 0, 0 ,0
#fig, ax = plt.subplots()
for delay, obsID in new_list_1:
    duration = duration_dict.get(obsID)
    truth = truth_dict.get(obsID)
    if duration is None:
        print(f"No duration found for obsID {obsID}")
        continue
    if truth is True:
        color = 'green'
        green = green + 1
    elif truth is False:
        color = 'red'
        red = red + 1
    else:
        color = 'grey'
        grey = grey + 1
sum = grey+red+green
print('this is red:', red, 'green:', green, 'grey:', grey, 'sum:',sum)
print('respective fractions', red/sum, green/sum, grey/sum)

   #ax.scatter(duration, delay, color=color)
#ax.set_xlabel('Duration of obsID (s)')
#ax.set_ylabel(r'Optimal delay $\tau$ with $C(\tau) > 3\sigma$ (s)')
#ax.set_ylim(0, 10000)
#plt.tight_layout()
#plt.show()


delays_1 = [delay for delay, obsID in new_list_1 if delay < 15000]
delays_2 = [delay for delay, obsID in new_list_2 if delay < 15000]


if plot_wanted_SPIE:
    #fig, ax = plt.subplots(figsize=(4.5, 6))
    extra = [d for d in delays_2 if d not in delays_1]  # or use set difference
    #ax.hist([delays_1, extra], bins = 10,stacked=True, color=['dodgerblue', '#84C7E6'],edgecolor='white', linewidth=0.5)
    #ax.hist(delays_2,  color="#84C7E6", edgecolor='white', linewidth=0.5)
    #ax.hist(delays_1,  color='dodgerblue', edgecolor='white', linewidth=0.5)
    #ax.set_xlabel(r'Optimal delays $\tau$ (s)', fontsize=12)
    #ax.set_ylabel('Number of correlated events', fontsize=12)
    #ax.tick_params(direction='in', length=6, width=1.5, labelsize=12, which='both')
    #ax.tick_params(length=3, which='minor')
    #ax.minorticks_on()
    #ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=6, prune='both'))
    #for spine in ax.spines.values():
     #   spine.set_linewidth(1.2)
    #plt.tight_layout()
    #plt.savefig('3sigma_histogram_newlist.png', dpi=150)
    #plt.show()


plot_wanted_final = False
if plot_wanted_final:
    #fig, ax = plt.subplots(figsize=(4.5, 6))
    fig, axs = plt.subplots(ncols = 1, nrows = 2)
    # ── main histogram ────────────────────────────────────────────────────────
    bins = np.histogram_bin_edges(np.concatenate([delays_1, delays_2]), bins=10)
    ax.hist([delays_1, extra], bins=bins, stacked=True,
            color=['dodgerblue', '#84C7E6'], edgecolor='white', linewidth=0.5)
    ax.set_xlabel(r'Optimal delays (s)', fontsize=12)
    ax.set_ylabel('Number of correlated events', fontsize=12)
    ax.tick_params(direction='in', length=6, width=1.5, labelsize=12, which='both')
    ax.tick_params(length=3, which='minor')
    ax.minorticks_on()
    ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=6, prune='both'))
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)

    # ── inset scatter ─────────────────────────────────────────────────────────
    # [x0, y0, width, height] in axes-relative coordinates
    ax_inset = ax.inset_axes([0.39, 0.56, 0.45*1.3, 0.30*1.3])

    for delay, obsID in new_list_1:
        duration = duration_dict.get(obsID)
        truth    = truth_dict.get(obsID)
        if duration is None:
            continue
        color = "forestgreen" if truth is True else ('firebrick' if truth is False else 'grey')
        ax_inset.scatter(duration, delay, color=color, s=18, linewidths=0)

    ax_inset.set_xlabel('ObsID length (s)', fontsize=11 )
    ax_inset.set_ylabel(r'$\tau$ (s)',  fontsize=11)
    ax_inset.set_ylim(0, 10000)
    ax_inset.tick_params(direction='in', length=3, width=1.0,
                         labelsize=10, which='both')
    ax_inset.tick_params(length=1.5, which='minor')
    ax_inset.minorticks_on()


    ax_inset.xaxis.set_major_locator(plt.MaxNLocator(nbins=4, prune='both'))
    #ax_inset.yaxis.set_major_locator(plt.MaxNLocator(nbins=3, prune='both'))

    for spine in ax_inset.spines.values():
        spine.set_linewidth(0.8)

    plt.tight_layout()
    plt.savefig('3sigma_histogram_newlist.png', dpi=150)
    plt.show()

if plot_wanted_final:
    fig, ax = plt.subplots(figsize=(4.5, 6))

    # ── main histogram ────────────────────────────────────────────────────────
    bins = np.histogram_bin_edges(np.concatenate([delays_1, delays_2]), bins=10)
    ax.hist(large_list_mod[1], bins=bins,
            color=['dodgerblue'], edgecolor='white', linewidth=0.5)
    ax.set_xlabel(r'Optimal delays (s)', fontsize=12)
    ax.set_ylabel('Number of correlated events', fontsize=12)
    ax.tick_params(direction='in', length=6, width=1.5, labelsize=12, which='both')
    ax.tick_params(length=3, which='minor')
    ax.minorticks_on()
    ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=6, prune='both'))
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)

    # ── inset scatter ─────────────────────────────────────────────────────────
    # [x0, y0, width, height] in axes-relative coordinates
    ax_inset = ax.inset_axes([0.39, 0.56, 0.45*1.3, 0.30*1.3])




    for truth, delay,duration in zip(reduced_truths, large_list[1],durations_of_enhancement_obsIDs):
        if truth == True:
            ax_inset.scatter(duration,delay, color = 'forestgreen')
        elif truth == False:
            ax_inset.scatter(duration, delay, color = 'firebrick')
        else:
            continue



    ax_inset.set_xlabel('ObsID length (s)', fontsize=11 )
    ax_inset.set_ylabel(r'$\tau$ (s)',  fontsize=11)
    ax_inset.set_ylim(0, 10000)
    ax_inset.tick_params(direction='in', length=3, width=1.0,
                         labelsize=10, which='both')
    ax_inset.tick_params(length=1.5, which='minor')
    ax_inset.minorticks_on()


    ax_inset.xaxis.set_major_locator(plt.MaxNLocator(nbins=4, prune='both'))
    #ax_inset.yaxis.set_major_locator(plt.MaxNLocator(nbins=3, prune='both'))

    for spine in ax_inset.spines.values():
        spine.set_linewidth(0.8)

    plt.tight_layout()
    plt.savefig('3sigma_histogram_newlist.png', dpi=150)
    plt.show()


plot_wanted_final_2 = True
if plot_wanted_final_2:
    fig, (ax2, ax1) = plt.subplots(ncols=2, nrows=1, figsize = (9,6))
    fs = 14
    fs_2 = 12
    # ── left histogram ────────────────────────────────────────────────────────
    bins = np.histogram_bin_edges(np.concatenate([delays_1, delays_2]), bins=20)
    ax1.hist([delays_1, extra], bins=bins, stacked=True,
             color=['dodgerblue', '#84C7E6'], edgecolor='white', linewidth=0.5)
    ax1.set_xlabel(r'Optimal delays (s)', fontsize=fs)
    ax1.set_ylabel('Number of correlated events', fontsize=fs)
    ax1.tick_params(direction='in', length=6, width=1.5, labelsize=fs, which='both')
    ax1.tick_params(length=3, which='minor')
    ax1.minorticks_on()
    ax1.set_xlim(xmin = 0)
    ax1.xaxis.set_major_locator(plt.MaxNLocator(nbins=6, prune='both'))
    for spine in ax1.spines.values():
        spine.set_linewidth(1.4)

    # left inset
    ax_inset1 = ax1.inset_axes([0.39, 0.46, 0.45*1.3, 0.40*1.3])
    for delay, obsID in new_list_1:
        duration = duration_dict.get(obsID)
        truth    = truth_dict.get(obsID)
        if duration is None:
            continue
        color = 'forestgreen' if truth is True else ('firebrick' if truth is False else 'grey')
        if color == 'forestgreen':
            ax_inset1.scatter(duration, delay, color=color, s=20, linewidths=0)
        if truth is False:
            ax_inset1.scatter(duration, delay, color=color, s=20,  marker = 'x')

    ax_inset1.set_xlabel('ObsID length (s)', fontsize=fs_2)
    ax_inset1.set_ylabel(r'$\tau$ (s)', fontsize=fs_2)
    ax_inset1.set_ylim(0, 15000)
    ax_inset1.tick_params(direction='in', length=3, width=1.2, labelsize=10, which='both')
    ax_inset1.tick_params(length=1.5, which='minor')
    ax_inset1.minorticks_on()
    ax_inset1.xaxis.set_major_locator(plt.MaxNLocator(nbins=4, prune='both'))
    for spine in ax_inset1.spines.values():
        spine.set_linewidth(1.2)

    # ── right histogram ───────────────────────────────────────────────────────
    ax2.hist(large_list_mod[1], bins=bins,
             color="#4E515B", edgecolor='white', linewidth=0.5, label = 'all delays')

    ax2.set_xlabel(r'Optimal delays (s)', fontsize=fs)
    ax2.set_ylabel('Number of correlated events', fontsize=fs)
    ax2.tick_params(direction='in', length=6, width=1.5, labelsize=fs, which='both')
    ax2.tick_params(length=3, which='minor')
    ax2.minorticks_on()
    ax2.set_xlim(xmin = 0)
    ax2.xaxis.set_major_locator(plt.MaxNLocator(nbins=6, prune='both'))
    for spine in ax2.spines.values():
        spine.set_linewidth(1.4)

    # right inset

    ax_inset2 = ax2.inset_axes([0.39, 0.46, 0.45*1.3, 0.40*1.3])
    for truth, delay, duration in zip(reduced_truths, large_list[1], durations_of_enhancement_obsIDs):
        if truth is True:
            ax_inset2.scatter(duration, delay, color='forestgreen', s = 12)
        elif truth is False:
            ax_inset2.scatter(duration, delay, color='firebrick', s = 16, marker = 'x')
    ax_inset2.set_xlabel('ObsID length (s)', fontsize=fs_2)
    ax_inset2.set_ylabel(r'$\tau$ (s)', fontsize=fs_2)
    ax_inset2.set_ylim(0, 15000)
    ax_inset2.tick_params(direction='in', length=3, width=1.2, labelsize=10, which='both')
    ax_inset2.tick_params(length=1.5, which='minor')
    ax_inset2.minorticks_on()
    ax_inset2.xaxis.set_major_locator(plt.MaxNLocator(nbins=4, prune='both'))
    for spine in ax_inset2.spines.values():
        spine.set_linewidth(1.2)


    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    legend_elements = [
        Patch(facecolor="#2E3139",  edgecolor='white', label='All delays'),
        Patch(facecolor='dodgerblue', edgecolor='white', label=r'$3\sigma$ delays'),
        Patch(facecolor='#84C7E6',    edgecolor='white', label = r'$3\sigma$ delays (including degeneracy)'),
        Patch(facecolor='royalblue',  edgecolor='white', label='Telescopes are on the:'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='forestgreen',
            markersize=7, label='Same side of the ecliptic'),
        Line2D([0], [0], marker='x', color='w', markerfacecolor='firebrick',
            markersize=7, label='Different sides of the ecliptic'),
    ]

    #fig.legend(handles=legend_elements, frameon=False, fontsize=10,
    #       loc='center left', bbox_to_anchor=(1.0, 0.5))

    fig.legend(handles=legend_elements, frameon=False, fontsize=10,
           loc='upper left', bbox_to_anchor=(0.92, 0.88),
           bbox_transform=fig.transFigure)

    #plt.savefig('_histograms.png', dpi=150, )
    #plt.show()


    #plt.tight_layout()
    plt.savefig('_histograms.png', dpi=150, bbox_inches='tight')
    plt.show()

