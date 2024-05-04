#!/usr/bin/python
import h5py
import matplotlib.pyplot as plt
import numpy as np
import math
import sys
from operator import itemgetter

number_of_pairs = 5

path_1 = "/PATH_OUTPOWERSCAN_ONE/POWERSCAN_ONE__int.h5" # INT
path_2 = "/PATH_OUTPOWERSCAN_TWO/POWERSCAN_TWO__int.h5"
# Distance threshold
epsilon = 1.0

# Open
fh_1 = h5py.File(path_1, 'r')
fh_2 = h5py.File(path_2, 'r')

n = 100
# FUNCTIONS
def dis (x, y):
	return [math.sqrt(xx**2 + yy**2) for xx, yy in zip(x, y)]

def column(matrix, i):
	return [row[i] for row in matrix]

def match(bitmask_1, bitmask_2):
	x1 = column(bitmask_1, 0)
	y1 = column(bitmask_1, 1)
	I1 = column(bitmask_1, 2)

	x2 = column(bitmask_2, 0)
	y2 = column(bitmask_2, 1)
	I2 = column(bitmask_2, 2)

	d1 = dis(x1, y1)
	d2 = dis(x2, y2)

	b1 = [[xx, yy, dd, II] for xx, yy, dd, II in zip(x1, y1, d1, I1)]
	b2 = [[xx, yy, dd, II] for xx, yy, dd, II in zip(x2, y2, d2, I2)]
	b1 = sorted(b1, key=itemgetter(2))
	b2 = sorted(b2, key=itemgetter(2))

	matched_pairs = []
	for i, p1 in enumerate(b1):
		for j, p2 in enumerate(b2):
			if(abs( p1[2]-p2[2] ) < epsilon):
				matched_pairs.append([p1, p2])
				del b2[j]
	return matched_pairs

def match2(matched_pairs_1, matched_pairs_2):
	matched_pairs_1_n = []
	matched_pairs_2_n = []
	for i, row1 in enumerate(matched_pairs_1):
		for j, row2 in enumerate(matched_pairs_2):
			if(abs(row1[0][2]-row2[0][2]) < epsilon):
				matched_pairs_1_n.append(row1)
				matched_pairs_2_n.append(row2)
				del matched_pairs_2[j]
	ratios_pairs = []
	for el1, el2 in zip(matched_pairs_1_n, matched_pairs_2_n):
		RAT1 = el1[0][3] / el1[1][3]
		RAT2 = el2[0][3] / el2[1][3]
		ratios_pairs.append([RAT1, RAT2])
	return ratios_pairs

# PROG
matched_pairs_1 = []
matched_pairs_2 = []

matched_pairs_1_n = []
matched_pairs_2_n = []
ratios_pairs = []
for angle  in range(41):
	print("ANGLE = {:d}".format(angle))
	for i in range (number_of_pairs):
		name_first =  "Reflection {:04d}".format(number_of_pairs * angle + i)
		name_second = "Reflection {:04d}".format(number_of_pairs * angle + i + 1)
		bitmask_1 = fh_1['Bitmasks Spots'][name_first]
		bitmask_2 = fh_1['Bitmasks Spots'][name_second]

		bitmask_3 = fh_2['Bitmasks Spots'][name_first]
		bitmask_4 = fh_2['Bitmasks Spots'][name_second]

		matched_pairs_1 = match(bitmask_1, bitmask_2)
		matched_pairs_2 = match(bitmask_3, bitmask_4)
		ratios_pairs = ratios_pairs + match2(matched_pairs_1, matched_pairs_2)


ratios_x = column(ratios_pairs, 0)
ratios_y = column(ratios_pairs, 1)


# Ploting
fonts_sizes = 10
fig = plt.figure()
ax = plt.gca()
im1 = ax.scatter(ratios_x, ratios_y, s=1)
ax.xaxis.set_tick_params(labelsize=fonts_sizes)
ax.yaxis.set_tick_params(labelsize=fonts_sizes)
plt.title('POWERSCAN_ONE-POWERSCAN_TWO', fontsize=fonts_sizes)
plt.xlim([0, 2])
plt.ylim([0, 2])
plt.gca().set_aspect('equal', adjustable='box')
plt.grid(True)
fig.savefig("/PATH_OUTcorr_plot_POWERSCAN_ONE-POWERSCAN_TWO.png", bbox_inches='tight', dpi=800)
#plt.show()

