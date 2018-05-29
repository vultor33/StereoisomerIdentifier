#READING AND ANALYZING RESULTS
import csv

k = 0
l = 0
with open('calculating.csv', 'r') as csvfile:
	spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')
	for row in spamreader:
		if len(row) > 3:
			if row[1] == 'Dimetal':
				stereoID1 = row[4].split('-')
				stereoID2 = row[8].split('-')
				if len(stereoID1) > 4 and len(stereoID2) > 4:
					k +=1
					if stereoID1[0] == stereoID2[0]:
						l+=1
						#print('1:  ',stereoID1[4], '  2:  ',stereoID2[4])
						#print('1:  ',stereoID1[3], '  2:  ',stereoID2[3])
					else:
						print(row)

print('total:  ',k)
print('same code: ',l)
