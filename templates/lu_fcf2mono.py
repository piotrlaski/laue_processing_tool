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
import optparse, sys

import numpy as np
import h5py

import laue_util

# -------------------------------------------------------------------------------------------------
usage = """
Program will read .fcf file and write hdf5 file with raw and processed mono data. \n
usage: %prog [options] file.fcf"""
parser = optparse.OptionParser(usage=usage)
parser.add_option("-o", "--output", type="str", dest="fo", default='mono.h5', 
                  help="(HDF5) output file name [default: %default]")
parser.add_option("-I", "--min_I", type="float", dest="I", default=0.0, 
                  help="minimum intensity [default: %default]")
parser.add_option("-r", "--max_h", type="float", dest="max_h", default=0.0, 
                  help="maximum |h|, 0.0 for no filtering  [default: %default]")
parser.add_option("-N", "--max_N", type="int", dest="max_N", default=0, 
                  help="maximum number of reflexions, 0 for all [default: %default]")
parser.add_option("-s", "--min_I_to_sigma", type="float", dest="min_I_to_sigma", default=0.0, 
                  help="minimum I/sigma, 0.0 for no filtering [default: %default]")
parser.add_option("-e", "--hkl_equivs", type="str", dest="hkl_equivs", default='[(h,k,l), (-h,-k,-l)]', 
                  help="laue symmetry [default: %default]. Should be enclosed in ''.")
parser.add_option("-c", "--cif", type="str", dest="cif", default='', 
                  help="optional: CIF file name. File content will be stored in the output file.")

(options, args) = parser.parse_args()
if len(args) != 1:
        parser.error("A single .fcf file name is a required argument.")
        
hkl_equivs = eval('lambda h,k,l: %s' % options.hkl_equivs)

# -------------------------------------------------------------------------------------------------

comment = """
FILE_CONTENTS:
 group      /
 dataset    /__snop_names__     
 group      /cell                 # unit cell parameters from fcf files
 dataset    /cell/A               # lattice vector a
 dataset    /cell/B               # lattice vector b
 dataset    /cell/C               # lattice vector c
 dataset    /cell/H               # reciprocal lattice vector a*
 dataset    /cell/K               # reciprocal lattice vector b*
 dataset    /cell/L               # reciprocal lattice vector c*
 dataset    /cell/U_HKL           # = array([H,K,L]).transpose()  -->  h = np.dot(hkl, U_HKL.transpose())
 dataset    /cell/V               # unit cell volume
 dataset    /cell/a               # |a|  (i.e. norm of lattice vector a)
 dataset    /cell/alpha           # alpha in radians
 dataset    /cell/alpha_deg       # alpha in degrees
 dataset    /cell/b
 dataset    /cell/beta
 dataset    /cell/beta_deg
 dataset    /cell/c
 dataset    /cell/c_star          # |c*|
 dataset    /cell/cos_alpha_star
 dataset    /cell/g_hkl           # = dot(U_HKL.transpose(), U_HKL)  -->  scalar_product(hkl1, hkl2) = np.dot(hkl1, np.dot(g_hkl, hkl2))
 dataset    /cell/gamma
 dataset    /cell/gamma_deg
 dataset    /fcf                  # full content of .fcf file
 group      /info
 dataset    /info/args
 dataset    /info/options
 group      /mono
 dataset    /mono/F_squared_calc
 dataset    /mono/F_squared_meas
 dataset    /mono/F_squared_sigma
 group      /mono/__attributes
 dataset    /mono/__attributes/MonoFilter_msg
 dataset    /mono/__attributes/cell
 dataset    /mono/__attributes/intensity            # choice of intensity for filtering and sorting (F_squared_meas or F_squared_calc)
 dataset    /mono/__class
 dataset    /mono/__keys
 dataset    /mono/h                                 # h vector in reciprocal space in cartesian coordinates
 dataset    /mono/h_len                             # |h|
 dataset    /mono/h_unit                            # h / |h|
 dataset    /mono/hkl                               # hkl : h = h * (a*) + k * (b*) + l * (c*),  h,k,l - integers
 dataset    /mono/idx_fcf                           # index in .fcf file
 dataset    /mono/idx_mono                          # index after filtering with MonoFilter i.e. in mono dataset                            
 group      /mono_fcf
 dataset    /mono_fcf/F_squared_calc
 dataset    /mono_fcf/F_squared_meas
 dataset    /mono_fcf/F_squared_sigma
 group      /mono_fcf/__attributes
 dataset    /mono_fcf/__attributes/cell
 dataset    /mono_fcf/__class
 dataset    /mono_fcf/__keys
 dataset    /mono_fcf/h
 dataset    /mono_fcf/h_len
 dataset    /mono_fcf/h_unit
 dataset    /mono_fcf/hkl
 dataset    /mono_fcf/idx_fcf
 group      /mono_x                                 # after MonoExpandFilter - contains symmetry equivalents
 dataset    /mono_x/F_squared_calc
 dataset    /mono_x/F_squared_meas
 dataset    /mono_x/F_squared_sigma
 group      /mono_x/__attributes
 dataset    /mono_x/__attributes/MonoFilter_msg
 dataset    /mono_x/__attributes/cell
 dataset    /mono_x/__attributes/intensity
 dataset    /mono_x/__class
 dataset    /mono_x/__keys
 dataset    /mono_x/h
 dataset    /mono_x/h_len
 dataset    /mono_x/h_unit
 dataset    /mono_x/h_unit_x
 dataset    /mono_x/h_x
 dataset    /mono_x/hkl                              # contains "merged" hkl-s i.e. as they are in mono
 dataset    /mono_x/hkl_x                            # contains symmetry generated hkl-s 
 dataset    /mono_x/idx_fcf
 dataset    /mono_x/idx_mono                         # index of corresponding entry in mono dataset
 dataset    /mono_x/idx_mono_x                       # index of entry in _this_ dataset
 dataset    /mono_x/in_fcf
 dataset    /mono_x/symop
 #
 # 
 # Rays are sets of hkl-s lying on a common half-line i.e.:
 #
 #  an hkl belongs to ray klm  <==>  exists an integer number gcd > 0, such that (h, k, l) = gcd * (k, l, m)  i.e.:
 #   h = gcd * k
 #   k = gcd * l 
 #   l = gcd * m
 #
 group      /rays_cc                                 # arrays with concatenated data for ray* to be indexed using *_BE
 group      /rays_cc/__attributes
 dataset    /rays_cc/__class
 dataset    /rays_cc/__keys
 dataset    /rays_cc/gcds                            # gratest common divisors:  (h, k, l) = gcd * (k, l, m) 
 dataset    /rays_cc/hkls                            # "merged" hkl-s belonging to a given ray
 dataset    /rays_cc/hklxs                           # symmetry generated hkl-s belonging to a given ray
 dataset    /rays_cc/idx_fcf                         # indexes of elements in fcf file contributing to a given ray
 dataset    /rays_cc/idx_mono                        # indexes of elements in mono dataset contributing to a given ray
 dataset    /rays_cc/idx_mono_x                      # indexes of elements in mono_x dataset contributing to a given ray
 group      /rays_fcf                                # rays belonging to space spanned by symmetry generated hkl-s
 group      /rays_fcf/__attributes
 dataset    /rays_fcf/__attributes/cell
 dataset    /rays_fcf/__class
 dataset    /rays_fcf/__keys
 dataset    /rays_fcf/gcds_BE
 dataset    /rays_fcf/h_len_at_max_F_squared_calc
 dataset    /rays_fcf/h_len_at_max_F_squared_meas
 dataset    /rays_fcf/hkls_BE
 dataset    /rays_fcf/hklxs_BE
 dataset    /rays_fcf/idx_ray_x
 dataset    /rays_fcf/idxs_fcf_BE
 dataset    /rays_fcf/idxs_mono_BE
 dataset    /rays_fcf/idxs_mono_x_BE
 dataset    /rays_fcf/k
 dataset    /rays_fcf/k_len
 dataset    /rays_fcf/k_unit
 dataset    /rays_fcf/klm
 dataset    /rays_fcf/max_F_squared_calc
 dataset    /rays_fcf/max_F_squared_meas
 dataset    /rays_fcf/min_h_len
 dataset    /rays_fcf/n
 group      /rays_half                  # rays with positive first nonzero coefficient i.e. k > 0 or k = 0 and l > 0 or ...
 group      /rays_half/__attributes
 dataset    /rays_half/__attributes/cell
 dataset    /rays_half/__class
 dataset    /rays_half/__keys
 dataset    /rays_half/gcds_BE
 dataset    /rays_half/h_len_at_max_F_squared_calc
 dataset    /rays_half/h_len_at_max_F_squared_meas
 dataset    /rays_half/hkls_BE
 dataset    /rays_half/hklxs_BE
 dataset    /rays_half/idx_ray_x
 dataset    /rays_half/idxs_fcf_BE
 dataset    /rays_half/idxs_mono_BE
 dataset    /rays_half/idxs_mono_x_BE
 dataset    /rays_half/k
 dataset    /rays_half/k_len
 dataset    /rays_half/k_unit
 dataset    /rays_half/klm
 dataset    /rays_half/max_F_squared_calc
 dataset    /rays_half/max_F_squared_meas
 dataset    /rays_half/min_h_len
 dataset    /rays_half/n
 group      /rays_x                  # all rays
 group      /rays_x/__attributes
 dataset    /rays_x/__attributes/cell
 dataset    /rays_x/__class
 dataset    /rays_x/__keys
 dataset    /rays_x/gcds_BE                        # [B:E] for /rays_cc/gcds
 dataset    /rays_x/h_len_at_max_F_squared_calc    # h_len of most intense (F_squared_calc) hkl contributing to the ray
 dataset    /rays_x/h_len_at_max_F_squared_meas    # h_len of most intense (F_squared_meas) hkl contributing to the ray
 dataset    /rays_x/hkls_BE
 dataset    /rays_x/hklxs_BE
 dataset    /rays_x/idx_ray_x
 dataset    /rays_x/idxs_fcf_BE
 dataset    /rays_x/idxs_mono_BE
 dataset    /rays_x/idxs_mono_x_BE
 dataset    /rays_x/k                              #  k (vector) = k * (a*) + l * (b*) + m * (c*)  (k,l,m - integers)
 dataset    /rays_x/k_len                          # |k|  
 dataset    /rays_x/k_unit                         # k / |K|
 dataset    /rays_x/klm                            # k,l,m integers
 dataset    /rays_x/max_F_squared_calc             # intensity of most intense (F_squared_calc) hkl contributing to the ray
 dataset    /rays_x/max_F_squared_meas             # intensity of most intense (F_squared_meas) hkl contributing to the ray
 dataset    /rays_x/min_h_len                      # minimum h_len over all hkl-s contributing to the ray
 dataset    /rays_x/n                              # number of hkls in the ray
"""

# -------------------------------------------------------------------------------------------------

mp = laue_util.mono_filter.MonoPipeline()
mp.mf.min_I = options.I
if options.max_h > 0.0:
    mp.mf.max_h_len = options.max_h
if options.max_N > 0.0:
    mp.mf.max_N = options.max_N
if options.min_I_to_sigma > 0.0:
    mp.mf.min_I_to_sigma = options.min_I_to_sigma
mp.mxf.hkl_equivs = hkl_equivs

fcf, txt = laue_util.read_fcf_hkls(args[0], return_raw_file=True)
mp.mf.mono_in = fcf
d = mp.get_snops() 
d['rays_cc'] = laue_util.snop.Snop(**mp.rf.cc_attributes)
print 'writing to %s' % options.fo
datasets = {'options': repr(options.__dict__), 'args': repr(args)}
f = laue_util.snop.write_snops(options.fo, **d)
laue_util.io.h5.mark_h5(f, '/home/piotr/Desktop/laueutil/laueutil-code/laue_util/bin/', datasets=datasets)
cg = f.create_group('cell')
for k,v in fcf.cell.iteritems():
    cg.create_dataset(k, data=v)
f.create_dataset('fcf', data=txt)
f.create_dataset('comment', data=comment)
f.create_dataset('hkl_equivs', data=options.hkl_equivs)
f.flush()
if options.cif:
    f.create_dataset('cif', data=file(options.cif).read())
f.close()
