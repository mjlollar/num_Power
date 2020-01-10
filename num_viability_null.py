#Compatible with Python 3.x
#Compatible with numpy 1.17.x

#Typical RAM 3 MB per 100 indivuals (for request)
#Mem request ~1000kb

import sys
import numpy as np
import scipy.stats as sp

file_name = sys.argv[1]
number_sequenced = int(sys.argv[2])

individuals = np.loadtxt(file_name, dtype=int, max_rows=number_sequenced)

#Define chromosomes(or arms) in 2D arrays
X = individuals[:,0:100]
twoL = individuals[:,100:200]
twoR = individuals[:,200:300]
threeL = individuals[:,300:400]
threeR = individuals[:,400:500]

#2D arrays with exclusion of focal chromosome(or arm) for comparison
noX = np.concatenate((twoL,twoR,threeL,threeR), axis=1)
noTwoL = np.concatenate((X,twoR,threeL,threeR), axis=1)
noTwoR = np.concatenate((X,twoL,threeL,threeR), axis=1)
noThreeL = np.concatenate((X,twoL,twoR,threeR), axis=1)
noThreeR = np.concatenate((X,twoL,twoR,threeL), axis=1)


# This function moves column-wise through array 1 and compares row values to  row 
# value of column in row 0 of array 2, then moves column-wise through array 1 and 
# compares to row values to those in row 1 of array 2, etc.

# One p-value is logged for each window comparison.

p_values = []

def chr_compare(array1, array2):
	for col_two in range(0, array2.shape[1]):
		for col_one in range(0, array1.shape[1]):
			count1 = []
			count2 = []
			count3 = []
			count4 = []
			i = 0
			for row in array1[:,col_one]: 
				try:
					if row == 0:
						if array2[i,col_two] == 2:
							count1.append(1)
						elif array2[i,col_two] != 2:
							count2.append(1)
						else:
							pass
					elif row == 1:
						if array2[i,col_two] == 2:
							count3.append(1)
						elif array2[i,col_two] != 2:
							count4.append(1)
						else:
							pass
					elif row == 2:
						if array2[i,col_two] == 2:
							count3.append(1)
						elif array2[i,col_two] != 2:
							count4.append(1)
						else:
							pass
					else:
						pass
					i += 1
				except IndexError:
					pass
					print("ERROR Caught")
			oddsratio, pvalue = sp.fisher_exact([[len(count1), len(count2)], [len(count3), len(count4)]])
			p_values.append(pvalue)
			count1.clear()
			count2.clear()
			count3.clear()
			count4.clear()

# Run Function for all combinations of chromosome arms (or X)

chr_compare(X,noX)
chr_compare(twoL, noTwoL)
chr_compare(twoR, noTwoR)
chr_compare(threeL, noThreeL)
chr_compare(threeR, noThreeR)

# Get the lowest p-value from all window comparisons and print it
# This will be the output of the file
min_pvalue = min(p_values)
print(min_pvalue)
