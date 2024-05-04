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
import sys, os

import numpy as np
import h5py

from enthought.traits.api import HasTraits, Instance, Range, Int, Float, Array, Property, \
                                 Bool, Button, Str, List, DelegatesTo, Array, Trait
from enthought.traits.ui.api import View, Item, Group, RangeEditor, TextEditor
from enthought.tvtk.tools import ivtk
from enthought.tvtk.api import tvtk

import laue_util
from laue_util.ui.sfp3d import _mk_fig, SingleFigurePanel

import jdVis.Data
import jdVis.gr as gr

import math

# -------------------------------------------------------------------------------------------------
Nb_arguments = len(sys.argv)-1
# Nb_arguments different to 3 or 4
if not(Nb_arguments == 4):
    print """
    # -------------------------------------------------------------------------------------------------

      usage: 
    
      ipython -wthread lu__hkl_assign.py expt.h5 mono.h5 file.mda infos.h5

    # -------------------------------------------------------------------------------------------------

    """
    raise SystemExit

file_expt = sys.argv[1]
file_mono = sys.argv[2]
file_mda = sys.argv[3]

lambda_cutoff = 0.01
lambda_N = 120

# -------------------------------------------------------------------------------------------------

mda = laue_util.read_mda(file_mda)                           # a dictionary of raw data from mda file
mda_f = laue_util.mda2lc.MdaToLC()
mda_f.lambda_cutoff = lambda_cutoff
mda_f.max_N = lambda_N
mda_f.mda_in = mda

f_expt, d_expt = laue_util.snop.read_snops(file_expt, do_eval=True) 
expt = d_expt['expt']

max_sin_t = max(expt.sin_t)
lambda_min = min(mda_f.lc.l)
max_h_len = 2.0 * max_sin_t / lambda_min
FrameIndices = list(f_expt['lauecollect/frame_number_global'])

f_mono, d_mono = laue_util.snop.read_snops(file_mono, do_eval=True) 
mono = d_mono['mono_x']
cell = mono.cell
mono_f = laue_util.mono_filter.MonoFilter()
mono_f.max_h_len = max_h_len
mono_f.mono_in = mono

# -------------------------------------------------------------------------------------------------

mono_lc = laue_util.mono_lc.MonoLC()
mono_lc.calculate_xyz = True
mono_lc.mono_f = mono_f
mono_lc.mda_f = mda_f
mda_f.edit_traits()

print 'mono_lc.N =', mono_lc.N 

if (Nb_arguments == 4):
   infosfile_name = sys.argv[4]
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
def are_same_ray(hkls):
    hkls = tuple(hkls)
    hkl1 = hkls[0]
    for hkl2 in hkls[1:]:
        gcd1, klm1 = laue_util.hkl.find_ray(hkl1)
        gcd2, klm2 = laue_util.hkl.find_ray(hkl2)
        if klm1 != klm2:
            return False
    return True

# -------------------------------------------------------------------------------------------------



class MAIN(SingleFigurePanel):
    
    # IN
    mono_f = Instance(laue_util.mono_filter.MonoFilter)
    mono_lc = Instance(laue_util.mono_lc.MonoLC)
    expt = Instance(laue_util.snop.Snop)
    #f_expt = h5py
    #f_mono = h5py 
    #file_mda = Str
    
    # params
    paic_cell_size = Float(0.004)
    cutoff = Float(0.002)
    #U_txt = Str
    #U = Array
    #Uinv = Array

    # output:
    N_not_assigned = Int
    N_single_hkl = Int
    N_single_ray = Int
    N_multiple_rays = Int

    # internal:
    U_type = Trait('Frame Rotation U', 'Global Rotation U', 'Frame Full U')
    update_button = Button('Calculate')
    show_3D = Bool(False)
    Sphere_radius = Float(0.001)
    file_name = Str('yyy.h5')
    save_button = Button('Save')
    GlobalRotationU = Array
    FrameRotationU = List
    FrameFullU = List
    GlobalRotationU_available = Bool
    FrameRotationU_available = Bool
    FrameFullU_available = Bool
    InvGlobalRotationU = Array
    InvFrameRotationU = List
    InvFrameFullU = List
    list_h_sint_rot = List
    v = Instance(ivtk.IVTK)

    # Selection of frame subset
    use_frame_selection = Bool
    selection_str = Str
    selection_list = List


    traits_view = View(Group(
                       Group(
                               Item('paic_cell_size', editor=TextEditor(evaluate=float, auto_set=False, enter_set=True)),
                               Item('cutoff', editor=TextEditor(evaluate=float, auto_set=False, enter_set=True)),
                               #Item('U_txt', editor=TextEditor(auto_set=False, enter_set=True)),
                               Item('use_frame_selection'),
                               Item('selection_str', enabled_when='(use_frame_selection)',editor=TextEditor(evaluate=str, auto_set=False, enter_set=True))
                           ), # main Group
                        Group(
                               Item('N_not_assigned', style='readonly'),
                               Item('N_single_hkl', style='readonly'),
                               Item('N_single_ray', style='readonly'),
                               Item('N_multiple_rays', style='readonly'),
#                               label = 'results',
                               show_border = True
                           ), # main Group   
#                        Item(name='msg', style='readonly', editor=TextEditor(multi_line=True), label='summary')
                        Item('U_type'),
                        Item('update_button', show_label=False),
                        Item('show_3D'),
                        Item('Sphere_radius', editor=TextEditor(evaluate=float, auto_set=False, enter_set=True)),
                        Item('file_name'),
                        Item('save_button', show_label=False)
                        ),
                        resizable=True,
                        title = 'hkl_assign_monolc'   
                        ) # View            


    def __init__(self, *args, **kw):
        super(MAIN, self).__init__(*args, **kw)

    def _update_button_fired(self):
        self.calculate()

    def _save_button_fired(self):

        self.expt['HKL'] = self.HKL
        self.expt['lambda'] = self.wavelength
        self.expt['flag_not_assigned'] = self.flag_not_assigned
        self.expt['flag_single_hkl'] = self.flag_single_hkl
        self.expt['flag_single_ray'] = self.flag_single_ray
        self.expt['flag_multiple_rays'] = self.flag_multiple_rays
        f = laue_util.snop.write_snops(self.file_name, expt=self.expt) 
        datasets = {'file_mda': os.path.abspath(file_mda), 
                    'args': repr(sys.argv),
                    'cutoff': self.cutoff}
        parents = [
            ('mono', os.path.basename(f_mono.filename), os.path.abspath(f_mono.filename), f_mono['info'].attrs['uuid']),
            ('expt', os.path.basename(f_expt.filename), os.path.abspath(f_expt.filename), f_expt['info'].attrs['uuid']),
        ]
        #laue_util.io.h5.mark_h5(f, __file__, datasets=datasets, parents=parents)
        laue_util.io.h5.mark_h5(f, 'lu__hkl_assign_monolc', datasets=datasets, parents=parents)
        print 'File written.'

    def calculate(self):
        if self.expt is None:
           raise ValueError, 'expt is not specified '
        if self.mono_f is None:
           raise ValueError, 'mono_f is not specified '
        if self.mono_lc is None:
           raise ValueError, 'mono_lc is not specified '

        if (self.U_type is None):
           print 'U type must be selected'
           return
        if (self.U_type == 'Global Rotation U') and not(self.GlobalRotationU_available == True):
           print 'Global rotation U matrix must be defined first'
           return
        if (self.U_type == 'Frame Rotation U') and not(self.FrameRotationU_available == True):
           print 'Frame rotation U matrices must be defined first'
           return
        if (self.U_type == 'Frame Full U') and not(self.FrameFullU_available == True):
           print 'Frame full U matrices must be defined first'
           return

        # Subset selection
        len_expt = len(self.expt)
        select_expt = np.ones(len_expt, dtype=bool)

        if self.use_frame_selection == True:
           self.process_selection()

        self.pcol = laue_util.PAiC.PAiC(cell_size=self.paic_cell_size)
        for i,x in enumerate(self.mono_lc.h_sint):
             self.pcol.add(i,x)
        #self.HKLs = []
        mono_indices = []
        list_h_sint_rot = []

        for i in xrange(len_expt):
            if (self.use_frame_selection == True):
               if not(self.expt.frame[i] in self.selection_list):
                  select_expt[i] = False
                  continue

	    h = np.dot(self.InvU_application(self.expt.frame[i]),self.expt.h_unit_g[i])
            h = h/math.sqrt(np.sum(h**2))
            h = h * self.expt.sin_t[i]
            list_h_sint_rot.append(h)
            shell = [idx for idx, y, d in self.pcol.shell(h, self.cutoff)]
            indices = set()
            for j in shell:
                mono_idx = j // self.mono_lc.N_lc
                #hkls.add(tuple(self.mono_lc.mono.hkl[mono_idx]))
                indices.add(mono_idx)
            #self.HKLs.append(hkls)
            mono_indices.append(indices)
        self.list_h_sint_rot = list_h_sint_rot
        #l = map(len, self.HKLs)
        l = map(len, mono_indices)
        print 'cutoff =', self.cutoff
        for i in set(l):
            #print i, len([s for s in self.HKLs if len(s) == i])
            print i, len([s for s in mono_indices if len(s) == i])
        self.HKL = np.zeros((len(self.expt),3), dtype='int')
        self.wavelength = np.zeros(len(self.expt), dtype='float')
        self.flag_not_assigned = np.zeros(len(self.expt), dtype='int')
        self.flag_single_hkl = np.zeros(len(self.expt), dtype='int')
        self.flag_single_ray = np.zeros(len(self.expt), dtype='int')
        self.flag_multiple_rays = np.zeros(len(self.expt), dtype='int')

        indices_expt = np.nonzero(select_expt)[0]
        #for i,s in enumerate(self.HKLs):
        for i,s in enumerate(mono_indices):
            hkls = [tuple(self.mono_lc.mono.hkl[elt]) for elt in list(s)]
            j = indices_expt[i]
            if len(s) == 0:
                self.flag_not_assigned[j] = 1
            elif len(s) == 1:
                self.flag_single_hkl[j] = 1
                self.HKL[j,:] = hkls[0]
                self.wavelength[j] = 2.0*self.expt.sin_t[j]/self.mono_lc.mono.h_len[list(s)[0]]
            elif are_same_ray(hkls):
                self.flag_single_ray[j] = 1
            else:
                self.flag_multiple_rays[j] = 1
        
        self.N_not_assigned = np.sum(self.flag_not_assigned)
        self.N_single_hkl = np.sum(self.flag_single_hkl)
        self.N_single_ray = np.sum(self.flag_single_ray)
        self.N_multiple_rays = np.sum(self.flag_multiple_rays)
        #for i,s in enumerate(self.HKLs):
        for i,s in enumerate(mono_indices):
            if self.flag_multiple_rays[i]: 
                l = [tuple(self.mono_lc.mono.hkl[j]) for j in s]
                print l
        self.update_viewer()

    def _show_3D_changed(self):
        self.viewer()

    def InvU_application(self, frame):
        if self.U_type == 'Global Rotation U':
           return self.InvGlobalRotationU
        if self.U_type == 'Frame Rotation U':
           return self.InvFrameRotationU[frame]
        if self.U_type == 'Frame Full U':
           return self.InvFrameFullU[frame]

# -------------------------------------------------------------------------------------------------
    def viewer(self):
        if (len(self.list_h_sint_rot) == 0):
           print 'hkl assigning must be performed first'
           return

        if (self.Sphere_radius <= 0):
           print 'Sphere radius must be positive'
           return

        if not(self.v is None): self.v.close()
        if self.show_3D:
           bonds = np.arange(self.mono_lc.N).reshape((self.mono_lc.N_mono, self.mono_lc.N_lc))
           p = jdVis.Data.PointsAndBonds(points=self.mono_lc.h_sint, bonds=bonds)
           p.add_atom_data(np.sqrt(self.mono_lc.I), 'sqrt_I')
           p.add_atom_data(self.mono_lc.l, 'l')

           ncolors = 256
           lut = laue_util.visualization.lut.hue_range_LUT(ncolors)

           pdgr = gr.SingleInputPolyDataGR(input=p.poly_data)
           pdgr.use_lookup_table_scalar_range = 1
           pdgr.lut.range = (min(self.mono_lc.l), max(self.mono_lc.l))
           p.poly_data.point_data.set_active_scalars('l') 

           #trf_expt = [h_rot * s for h_rot, s in zip(self.list_h_sint_rot, self.expt.sin_t)]
           trf_expt = self.list_h_sint_rot
           pe = jdVis.Data.PointsAndBonds(points=trf_expt)
           spheres1 = gr.Spheres(input=pe.poly_data, radius=self.Sphere_radius)
           spheres1.theta_resolution = 4
           spheres1.phi_resolution = 4

           self.v = ivtk.IVTK(size=(800,800))
           self.v.open()
           self.v.scene._renwin.stereo_type = 'anaglyph'
           self.v.scene._renwin.anaglyph_color_mask = 4,2
           # here are lines (RK)
           #self.v.scene.add_actor(pdgr.actor)
           self.v.scene.add_actor(spheres1.actor)
	   self.spheres = spheres1

    def update_viewer(self):
        if self.show_3D:
           if (self.Sphere_radius <= 0):
              print 'Sphere radius must be positive'
              return
	   #trf_expt = [h_rot * s for h_rot, s in zip(self.list_h_sint_rot, self.expt.sin_t)]
           trf_expt = self.list_h_sint_rot
           pe = jdVis.Data.PointsAndBonds(points=trf_expt)
           spheres1 = gr.Spheres(input=pe.poly_data, radius=self.Sphere_radius)
           spheres1.theta_resolution = 4
           spheres1.phi_resolution = 4

           self.v.scene.remove_actor(self.spheres.actor)
           self.v.scene.add_actor(spheres1.actor)
	   self.spheres = spheres1

    def _Sphere_radius_changed(self):
	   self.update_viewer()

    def process_selection(self):
        selection_str = self.selection_str.replace(' ','').replace(';',',') 
        selection_list = []
        if (len(self.selection_str) > 0):
           list_words = selection_str.split(',')
           for word in list_words:
              if ('-' in word):
                 list_indices = word.split('-')
                 if (len(list_indices) == 2):
                    if (int(list_indices[0]) <= int(list_indices[1])):
                       selection_list.extend(range(int(list_indices[0]),int(list_indices[1])+1))
                 else:
                    print 'Error in subset definition: ', word
              else:
                 if not(word == ''):
                    selection_list.append(int(word))

        # Remove redundances
        selection_list = sorted(list(set(selection_list) & set(FrameIndices)))
        if (len(selection_list)==0):
           selection_list = FrameIndices[:]
        self.selection_list = selection_list

        # Rewrite the selection string
        new_selection_str = ''
        len_selection_list = len(selection_list)
        if (len_selection_list==1):
           new_selection_str = str(selection_list[0])
        else:
           for i, indice in enumerate(selection_list):
              if (i == 0):
                 indice_begin = selection_list[0]
                 previous_indice = indice_begin
                 continue
              go_on = False
              if (indice == previous_indice+1):
                 go_on = True
              if (go_on == True):
                 previous_indice = indice
              else:
                 str_size = len(new_selection_str)
                 if not(str_size == 0):
                    new_selection_str = new_selection_str + ','
                 if (indice_begin == previous_indice):
                    new_selection_str = new_selection_str + str(indice_begin)
                 else:
                    new_selection_str = new_selection_str + str(indice_begin) + '-' + str(previous_indice)
                 previous_indice = indice
                 indice_begin = indice

              if (len_selection_list == i+1):
                 str_size = len(new_selection_str)
                 if not(str_size == 0):
                    new_selection_str = new_selection_str + ','
                 if (indice_begin == previous_indice):
                    new_selection_str = new_selection_str + str(indice_begin)
                 else:
                    new_selection_str = new_selection_str + str(indice_begin) + '-' + str(previous_indice)
                 
        self.selection_str = new_selection_str
           
m = MAIN()
m.mono_f = mono_f
m.mono_lc = mono_lc
m.expt = expt
m.use_frame_selection = False
m.selection_str = ''

infosfile = h5py.File(infosfile_name, 'r')

# Check if initial matrix evaluable in infos_file
if ('Global_Rotation_Matrix' in infosfile.listnames()):
   m.GlobalRotationU_available = True
   m.GlobalRotationU = np.array(infosfile['Global_Rotation_Matrix/Indirect_Rotation_Matrix'])
   m.InvGlobalRotationU = np.linalg.inv(m.GlobalRotationU)
else:
   m.GlobalRotationU_available = False

# Check if frame rotation matrices evaluable in infos file
if ('Frame_Rotation_Matrices' in infosfile.listnames()):
   m.FrameRotationU_available = True
   AllFrameRotationU = np.array(infosfile['Frame_Rotation_Matrices/Indirect_Rotation_Matrices'])
   m.FrameRotationU = [AllFrameRotationU[k,:].reshape((3,3), order='C') for k in xrange(AllFrameRotationU.shape[0])]
   m.InvFrameRotationU = [np.linalg.inv(FrameMatrix) for FrameMatrix in m.FrameRotationU]
else:
   m.FrameRotationU_available = False

# Check if frame transformation  matrices evaluable in infos file
if ('Frame_Transformation_Matrices' in infosfile.listnames()):
   m.FrameFullU_available = True
   AllFrameFullU = np.array(infosfile['Frame_Transformation_Matrices/Indirect_Transformation_Matrices'])
   m.FrameFullU = [AllFrameFullU[k,:].reshape((3,3), order='C') for k in xrange(AllFrameFullU.shape[0])]
   m.InvFrameFullU = [np.linalg.inv(FrameMatrix) for FrameMatrix in m.FrameFullU]
else:
   m.FrameFullU_available = False

infosfile.close()

m.edit_traits()

import os
import time

time.sleep(3)
m.calculate()
time.sleep(3)
try:
   m._save_button_fired()
except:
   pass
print 'End'

exit()
