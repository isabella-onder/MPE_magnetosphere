from datetime import datetime, timedelta

global starting_point 
starting_point = datetime(1999,12,31,21,00,00)


with open("long_revno", "r") as fp:
    data = fp.readlines()
    big_list = []
    for line in data[3589:]:
        parts = line.split()
        timestamp_xmm = parts[3]
        revs = int(parts[6])
        timestamp_xmm_ETC = datetime.strptime(timestamp_xmm, "%Y-%m-%dT%H:%M:%SZ")
        big_list.append((timestamp_xmm_ETC, revs))

#print(big_list)

def eRODAY_to_time_ETC(eroday):
    hours_passed = eroday * 4
    #starting_point = datetime(1999,12,31,21,00,00)
    time_eROSITA_ETC = starting_point + timedelta(hours=(hours_passed - 4))
    return time_eROSITA_ETC


#print(eRODAY_to_time_ETC(1))


def time_ETC_to_eRODAY(time_eROSITA_ETC):
    time_eROSITA_ETC = datetime.fromisoformat(time_eROSITA_ETC)
    print(time_eROSITA_ETC)
    hours_passed = (time_eROSITA_ETC - starting_point).total_seconds() / 3600 + 4
    eroday = hours_passed / 4
    return eroday

#print(time_ETC_to_eRODAY("2019-09-03T01:06:00"))



def revs_to_eROday(rev):
    timestamp = None 
    for i in range (len(big_list)):
        if rev == big_list[i][1]:
            timestamp = big_list[i][0]
            break
    if timestamp == None:
        raise ValueError("The timestamp was not found")

    time_difference = (timestamp - starting_point).total_seconds()/3600
    #print("these are the ones that do work", type(timestamp - starting_point), type(time_difference))
    time_difference_eROday = time_difference/4 +1
    return time_difference_eROday    

#print(revs_to_eROday())

def eROday_to_revs(eROday):
    timestamp = eRODAY_to_time_ETC(eROday)
    timestamp_seconds = (timestamp - datetime(2000,1,1,0,0,0)).total_seconds()
    seconds_column = [(row[0]-datetime(2000,1,1,0,0)).total_seconds() for row in big_list]

    i = 0
    while timestamp_seconds > seconds_column[i]:
        i +=1
    revs = big_list[i-1][1]
    return revs

print(eROday_to_revs(43205), eROday_to_revs(43206), eROday_to_revs(43215))



    

'''
Plans for tomorrow:
- revs to eroday:
       once given a rev, need to search through column and find corresponding timestamp. Convert it to datettime, find
       the delta t with the moscou 0 time for erosita in hours, and then find the corresponding eroday
- eroday to revs:
        once we have the timestamp for eroday, need to seach in the list of timestamps we made from the xmm ascii file - go 
        direct for date and the +/- to find correct day and hence rev to find closest timestamp
        to do so can just convert everything to seconds and then encadrer our value to find the corresponding rev
'''

'''
Alternative way to write it (whereas "with" has the advantage of automatically closing the file (fp.close()) when an error occurs)

fp = open("long_revno", "r") 
fp.read()
...
fp.close()
'''




