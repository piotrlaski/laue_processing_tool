# ---------------------------------------------------------------------------- #
#                                                                              #
#            LaueUtil  -  http://sourceforge.net/projects/laueutil/            #
#                                                                              #
# Developed at Philip Coppens Lab                                              #
#                                                                              #
# written by:                                                                  #
#                                                                              #
#     Jaroslaw A. Kalinowski  (2010-2011)  http://jak.kalinowscy.eu            #
#     Bertrand G.M. Fournier  (2011-2013)  betrandf@buffalo.edu                #
#                                                                              #
# ---------------------------------------------------------------------------- #

# The contents of this file are subject to the University at Buffalo Public 
# License Version 1.0 (the "License"); you may not use this file except in
# compliance with the License. You may obtain a copy of the License at
# http://sourceforge.net/projects/laueutil/files/License/ .
#
# Software distributed under the License is distributed on an "AS IS" basis,
# WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License 
# for the specific language governing rights and limitations under the License.
#
# The Original Code is LaueUtil software suite.
#
# The Initial Developer of the Original Code is Research Foundation of State
# University of New York, on behalf of University at Buffalo.
#
# Portions created by the Initial Developer are Copyright (C) 2007 Research
# Foundation of State University of New York, on behalf of University at 
# Buffalo. All Rights Reserved.
#
# Contributor(s): ______________________________________.
#
# Alternatively, the contents of this file may be used under the terms of either
# the GNU General Public License Version 2 (the "GPL"), or the GNU Lesser
# General Public License Version 2.1 (the "LGPL"), in which case the provisions
# of the GPL or the LGPL are applicable instead of those above. If you wish to
# allow use of your version of this file only under the terms of either the GPL
# or the LGPL, and not to allow others to use your version of this file under 
# the terms of the UBPL, indicate your decision by deleting the provisions above
# and replace them with the notice and other provisions required by the GPL or
# the LGPL. If you do not delete the provisions above, a recipient may use your
# version of this file under the terms of any one of the UBPL, the GPL or the 
# LGPL.

# ------------------------------------------------------------------------------
import optparse, os, sys, pprint, time
import numpy as np

import h5py
import laue_util

pp = pprint.PrettyPrinter(indent=4, width=120)

# -------------------------------------------------------------------------------------------------
usage = "usage: %prog [options] indexing.h5"
parser = optparse.OptionParser(usage=usage)
parser.add_option("-r", "--ratios_output", type="str", dest="ratios_output", default='LASER__ratios.h5',
                  help="output hdf5 file name for all ratios")
parser.add_option("-R", "--ratios_txt_output", type="str", dest="ratios_txt_output", default='LASER__ratios.txt',
                  help="output txt file name for all ratios")
parser.add_option("-o", "--output", type="str", dest="output", default='LASER__ratios.hkl',
                  help="output file name (sortav compatible)")
parser.add_option("-s", "--max_std", type="float", dest="max_std", default=0.5,
                  help="block filter: maximum sample standard deviation [default: %default]")
parser.add_option("-n", "--min_n", type="int", dest="min_n", default=1,
                  help="block filter: minimum number of pairs in block [default: %default]")
parser.add_option("-i", type="float", dest="i", default=0.0,
                  help="pair filter: I1 >= i and I2 >= i [default: %default]")
parser.add_option("-I", type="float", dest="minI", default=0.0,
                  help="minimum <Ioff> in block [default: %default]")
parser.add_option("-S", "--strategy", type="str", dest="strategy", default="ONOFF",
                  help="strategy used during laue expriments ONOFF or OFFON [default: %default]")
parser.add_option("-c", type="float", dest="c", default=0.0,
                  help="pair filter: I1 >= c*sig(I2) and I2 >= c*sig(I2) [default: %default]")
parser.add_option("--ddof", type="int", dest="ddof", default=1,
                  help="see nupy.std documentation [default: %default]")
parser.add_option("-P", "--partition", type="str", dest="partition", default="none",
                  help="partition strategy of eta values: 'none','angle','delay' [default: %default]")
parser.add_option("-l", "--lambdamin", type="float", dest="lambdamin", default=None,
                  help="Lambda min [default: %default]")
parser.add_option("-L", "--lambdamax", type="float", dest="lambdamax", default=None,
                  help="Lambda max [default: %default]")

(options, args) = parser.parse_args()
if len(args) == 0: parser.error("hkl indexing h5 file must be provided.")
tm = time.strftime("%Y-%m-%d-%H:%M:%S-%Z")
if len(options.output):        output = options.output
else:                          output = 'lu_indexing2ratios__%s__filtered.dat' % \
                                           (tm,)
                                        #   (tm, options.max_std, options.min_n, options.a, options.i, options.ddof)
if len(options.ratios_output): ratios_output = options.ratios_output
else:                          ratios_output = 'lu_indexing2ratios__%s__all.h5' % (tm,)
if len(options.ratios_txt_output): ratios_txt_output = options.ratios_txt_output
else:                          ratios_txt_output = 'lu_indexing2ratios__%s__all.dat' % (tm,)

if len(options.partition):
   partition = options.partition
if not(partition in ['none','delay','angle']):
   partition = 'none'
available_delays = False
lambdamin = options.lambdamin 
lambdamax = options.lambdamax
# -------------------------------------------------------------------------------------------------
def indexing2ratios(indexing_h5, max_std, min_n, min_i, min_I, sigmac, strategy, available_delays, ddof):

    # Selecting HKL assigned only
    pairs_of_frames = np.reshape(np.array(list(set(indexing_h5['expt/frame']))),(-1, 2),order='C') #/ON OFF frames/
    nb_pairs = pairs_of_frames.shape[0]
    print 'nb of pairs ', nb_pairs 
    if (strategy=='OFFON'): pairs_of_frames = pairs_of_frames[:, [1,0]]

    # Extract information
    D_R = {}   # key = hkl, angle, series  --> ratios 
    D_Ioff = {}   # key = hkl, angle, series  --> ratios
    D_Ion = {}   # key = hkl, angle, series  --> ratios
    list_of_ratios = []
    list_delays = []
    list_angles = []
    for pair_num in xrange(nb_pairs):
        # Searching common hkl between frames ON and OFF
        frames_num = pairs_of_frames[pair_num,:]
        mask_frameON = np.logical_and(np.array(indexing_h5['expt/flag_single_hkl'])==1, np.array(indexing_h5['expt/frame']) == frames_num[0])
        mask_frameOFF = np.logical_and(np.array(indexing_h5['expt/flag_single_hkl'])==1, np.array(indexing_h5['expt/frame']) == frames_num[1])

        # set of hkl
        if (np.shape(np.nonzero(mask_frameON)[0])[0]==0) or (np.shape(np.nonzero(mask_frameOFF)[0])[0]==0):
           continue
        hkl_frameON = indexing_h5['expt/HKL'][mask_frameON,:]
        hkl_frameOFF = indexing_h5['expt/HKL'][mask_frameOFF,:]
        List_hkl_frameON = map(tuple, hkl_frameON)
        List_hkl_frameOFF = map(tuple, hkl_frameOFF)
        d2 = dict(zip(List_hkl_frameON, np.nonzero(mask_frameON)[0]))
        d1 = dict(zip(List_hkl_frameOFF, np.nonzero(mask_frameOFF)[0]))

        Intersection_set = set(List_hkl_frameON) & set(List_hkl_frameOFF)
        for hkl in Intersection_set:
            i1 = d1[hkl]
            i2 = d2[hkl]
            hkltxt = '%+03d%+03d%+03d' % hkl
            hkltxt = hkltxt.replace('-','m').replace('+','p')
            h,k,l = hkl 
            I1 = indexing_h5['expt/I'][i1]
            I2 = indexing_h5['expt/I'][i2]
            angle = indexing_h5['expt/phi'][i1]
            Lambda = indexing_h5['expt/lambda'][i1]
            sin_t = indexing_h5['expt/sin_t'][i1]
            if not(lambdamax is None):
               if (Lambda > lambdamax): continue
            if not(lambdamin is None):
               if (Lambda < lambdamin): continue

            if not(angle in list_angles):
               list_angles.append(angle)
            if (available_delays == True):
               delay = indexing_h5['expt/delay'][i2]
               if not(delay in list_delays):
                  list_delays.append(delay)
            else:
               delay = 0.0

            # adding of peak coordinates
            XY1 = indexing_h5['expt/xy_raw'][i1]
            XY2 = indexing_h5['expt/xy_raw'][i2]
            if I1 > min_i and I2 > min_i:    
                ratio_data = \
                      (
                       frames_num[1], frames_num[0], 
                       hkltxt, 
                       h,k,l,
                       angle,delay,sin_t,Lambda, 
                       I1, I2,
                       I2 / I1,
                       XY1, XY2)
                list_of_ratios.append(ratio_data)
                key = hkl, angle, delay
                if not D_R.has_key(key):
                    D_R[key] = []
                if not D_Ioff.has_key(key):
                    D_Ioff[key] = []           
                if not D_Ion.has_key(key):
                    D_Ion[key] = []           
                D_R[key].append(I2 / I1)
                D_Ioff[key].append(I1)
                D_Ion[key].append(I2)
    KEYS = D_R.keys()
    KEYS.sort()
    
    D2 = {}
    for key in KEYS:
        hkl, angle, delay = key
        a = np.array(D_R[key])
        h,k,l = hkl 
        avg_a = np.average(a)
        std_a = np.std(a, ddof=ddof)
        avgI = np.average(D_Ioff[key])
        avgIoff = np.average(D_Ioff[key])
        avgIon = np.average(D_Ion[key])
        stdIoff = np.std(D_Ioff[key],ddof=ddof)
        stdIon = np.std(D_Ion[key], ddof=ddof)
        if avgI < min_I:
            continue
        if a.shape[0] >= min_n and std_a <= max_std and avgIon >= sigmac*stdIon and avgIoff >= sigmac*stdIoff:
            data = (h,k,l, avg_a, std_a)
            D2[key] = data
    
    return list_of_ratios, D2, list_angles, list_delays

# -------------------------------------------------------------------------------------------------

indexing_h5 = h5py.File(args[0], 'r')
if (not options.strategy in ['ONOFF','OFFON']):
   options.strategy = 'ONOFF'

if ('delay' in indexing_h5['expt'].keys()):
   available_delays = True

list_of_ratios, D2, list_angles, list_delays = indexing2ratios(indexing_h5,
                                    max_std = options.max_std, 
                                    min_n = options.min_n, 
                                    min_i = options.i,
                                    min_I = options.minI,
                                    sigmac = options.c,
                                    strategy = options.strategy,
                                    available_delays = available_delays,
                                    ddof = options.ddof)

f = h5py.File(ratios_output,'w')
f_txt = file(ratios_txt_output, 'w')

f_pairs = f.create_group("I_pairs")
total_nb_pairs = len(list_of_ratios)
I1_array = np.zeros(total_nb_pairs)
I2_array = np.zeros(total_nb_pairs)
R_array = np.zeros(total_nb_pairs)
angle_array = np.zeros(total_nb_pairs)
delay_array = np.zeros(total_nb_pairs)
lambda_array = np.zeros(total_nb_pairs)
sint_array = np.zeros(total_nb_pairs)
frame1_array = np.zeros(total_nb_pairs,dtype=np.int16)
frame2_array = np.zeros(total_nb_pairs,dtype=np.int16)
XY1_array = np.zeros((total_nb_pairs,2),dtype=np.int16)
XY2_array = np.zeros((total_nb_pairs,2),dtype=np.int16)
hkl_array = np.zeros((total_nb_pairs,3),dtype=np.int16)
index = 0
if (available_delays == True):
   txt = ' fr1 fr2    h    k    l     phi   delay   sin_t  lambda      int1      int2     ratio X_fr1 Y_fr1 X_fr2 Y_fr2\n'
else:
   txt = ' fr1 fr2    h    k    l     phi   sin_t  lambda      int1      int2     ratio X_fr1 Y_fr1 X_fr2 Y_fr2\n'
f_txt.write(txt)      
for ratio_data in list_of_ratios:
    # adding of peak coordinates
    i,j, hkltxt, h,k,l, angle, delay, sin_t, Lambda, I1, I2, R, XY1, XY2 = ratio_data
    I1_array[index] = I1
    I2_array[index] = I2
    R_array[index] = R
    angle_array[index] = angle
    delay_array[index] = delay
    lambda_array[index] = Lambda
    sint_array[index] = sin_t
    frame1_array[index] = i
    frame2_array[index] = j
    XY1_array[index,:] = np.transpose(XY1)
    XY2_array[index,:] = np.transpose(XY2)
    hkl_array[index,:] = np.array([h,k,l])
    index += 1
    # txt output file
    if (available_delays == True):
       txt = '%4d%4d%5d%5d%5d%8.3f%8.3f%8.5f%8.5f%10.3f%10.3f%10.3f%6d%6d%6d%6d\n' % (i,j,h,k,l,angle,delay,sin_t,Lambda,I1,I2,R,XY1[0],XY1[1],XY2[0],XY2[1])
    else:
       txt = '%4d%4d%5d%5d%5d%8.3f%8.5f%8.5f%10.3f%10.3f%10.3f%6d%6d%6d%6d\n' % (i,j,h,k,l,angle,sin_t,Lambda,I1,I2,R,XY1[0],XY1[1],XY2[0],XY2[1])
    f_txt.write(txt)      
    
f_pairs.create_dataset('hkl', data=hkl_array)
f_pairs.create_dataset('I1', data=I1_array)
f_pairs.create_dataset('I2', data=I2_array)
f_pairs.create_dataset('R', data=R_array)
f_pairs.create_dataset('angle', data=angle_array)
f_pairs.create_dataset('lambda', data=lambda_array)
f_pairs.create_dataset('sint', data=sint_array)
if (available_delays == True): f_pairs.create_dataset('delay', data=delay_array)
f_pairs.create_dataset('frame1', data=frame1_array)
f_pairs.create_dataset('frame2', data=frame2_array)
f_pairs.create_dataset('XY1', data=XY1_array)
f_pairs.create_dataset('XY2', data=XY2_array)
f.close()
f_txt.close()

KEYS = D2.keys()
KEYS.sort()
f = file(output, 'w')

for key in KEYS: 
    data = D2[key]
    txt = '%5d%5d%5d%15.7f%15.7f' % data
    if (partition=='none'):
       txt2 = '    1'
    if (partition=='delay'):
       index = list_delays.index(key[2])+1
       txt2 = '   %2d' % index
    if (partition=='angle'):
       index = list_angles.index(key[1])+1
       txt2 = '   %2d' % index
    txt = txt+txt2+'\n'
    f.write(txt)      
f.close()
print 'End'
time.sleep(1)
exit()
