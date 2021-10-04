import sqlite3
import os
home_path = os.path.expanduser('~/')
path_67p = home_path+'Documents/year1/shape_modelling/67p/'

conn = sqlite3.connect('baseline_v1.3_10yrs.db')
cur = conn.cursor()

print('\nColumns in summaryAllProps table:')
data=cur.execute('''SELECT * FROM SummaryAllProps''')
for column in data.description:
    print(column[0])


obs_ids = []
file = open(path_67p+'LSST_out_LC_67p.txt')
lines = file.readlines()
for line in lines[25:]:
    fata = line.split()
    obs_ids.append(str(fata[1]))
print(len(obs_ids))
# Display data

#string+=observationId+","
searchString="select fiveSigmaDepth, observationId from summaryAllProps where observationId in ("

for item in obs_ids:
    searchString+=item
    if (obs_ids.index(item)<len(obs_ids)-1):
        searchString+=","
    else:
        searchString+=")"

print('\nData in summaryAllProps table:')
data = cur.execute(searchString)
#data=cur.execute('''SELECT fiveSigmaDepth, observationId FROM SummaryAllProps where observationId=''')
for row in data:
    print(row)



conn.close()
