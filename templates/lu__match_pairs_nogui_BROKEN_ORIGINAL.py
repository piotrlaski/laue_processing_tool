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
import imp, sys

import numpy as np
import matplotlib.pyplot as plt

from enthought.traits.api import HasTraits, Instance, Int, Float, Button, Str, Bool, List, Array, \
                                 Range, Property, DelegatesTo, Trait
from enthought.traits.ui.api import View, Item, Group, TextEditor, RangeEditor
from enthought.tvtk.tools import ivtk
from enthought.tvtk.api import tvtk

import transformations

import laue_util
import h5py, copy

import jdVis.gr as gr
from jdVis.Data import PointsAndBonds

# -------------------------------------------------------------------------------------------------

if len(sys.argv) != 3:
    msg = '''
***************************************************************************************************
    
    Usage:
    
    ipython -wthread lu__match_pairs_clustering.py file_expt.h5 file_mono.h5
    
***************************************************************************************************
    '''
    sys.exit(msg)

# -------------------------------------------------------------------------------------------------
# Expt data
fe, de = laue_util.snop.read_snops(sys.argv[1])
expt = de['expt']
# Monochromatic data
fm, dm = laue_util.snop.read_snops(sys.argv[2])
rays = dm['rays_x']
mono = dm['mono_x']

Totalnbframes = fe['lauecollect/frame_number_global'].shape[0]
FrameIndices = fe['lauecollect/frame_number_global']

# 
# -------------------------------------------------------------------------------------------------

expt_f = laue_util.expt_filter.ExptFilter()
expt_f.block_size = fe['lauecollect'].attrs.get('block_size', 0)
if expt_f.block_size:
    expt_f.position_in_block = 0
expt_f.expt_in = expt

ecf = laue_util.expt_filter.ExptClusteringFilter(use_idx_expt=False, min_I=50000)
ecf.expt_f = expt_f

ecm = laue_util.ui.expt_clusters.ExptClustersManager()
ecm.ecf = ecf

# -------------------------------------------------------------------------------------------------

expt_f.edit_traits()
c1_expt_ui = ecf.c1.edit_traits()
c1_expt_ui.title = 'Cluster Expt'
ecm.edit_traits()
ecv = laue_util.ui.expt_clusters.ExptClustersViewer()
ecv.ecm = ecm
ecv.edit_traits()
ecv.update()
            
# -------------------------------------------------------------------------------------------------

rf = laue_util.mono_filter.RayFilter(max_nb_selected=100)
rf.rays_in = rays

rf.edit_traits()

# -------------------------------------------------------------------------------------------------

class Picker(HasTraits):

    def __init__(self, scene, spheres, handler):
        self.scene = scene
        self.spheres = spheres
        self.handler = handler
        self.cellpicker = tvtk.CellPicker()
        self.cellpicker.tolerance = 0.0
    
    def pick(self, x, y):
        self.cellpicker.pick(float(x), float(y), 0.0, self.scene.renderer)
        if self.cellpicker.cell_id > -1:
            pd = self.spheres.tv_glyph.get_output()
            if self.cellpicker.data_set is pd:
                point_id = pd.get_cell(self.cellpicker.cell_id).point_ids[0]
                index = int(pd.point_data.get_array("InputPointIds")[point_id])
                self.handler(index)

def format_inp(crystal, U):
    txt = '''Input
   Crystal    %s  
   Matrix     %s
'''
    crystal = '%f %f %f   %f %f %f' % crystal
    mtx = '%f %f %f  ' % tuple(U[1])
    mtx += '%f %f %f  ' % tuple(U[2])
    mtx += '%f %f %f' % tuple(U[0])
    return txt  % (crystal, mtx)

class Main(HasTraits):
    
    angular_cutoff = Float(0.5)
    #txt = Str('')
    show_3D = Bool(False)
    calculate_button = Button('Calculate')
    cluster_button = Button('Cluster')
    orient_matrix_selection = Int
    infos_filename = 'xxx.h5'
    save_button = Button('Save exptfile')
    UsAreAvailable = False 
    c1 = Instance(laue_util.cluster_points.Cluster1)
    v = None

    main_group = Group(
        Item('angular_cutoff', editor=TextEditor(evaluate=float, auto_set=False, enter_set=True)),
        #Item('txt', editor=TextEditor(evaluate=float, auto_set=False, enter_set=True)),
        Item('show_3D'),
        Item('calculate_button', show_label=False),
        Item('cluster_button', show_label=False),        
        Item('orient_matrix_selection', editor=TextEditor(evaluate=int, auto_set=False, enter_set=True)),
        Item('infos_filename'),
        Item('save_button', show_label=False),
        label = 'Main',
        show_border = False
    )
    traits_view = View(main_group, resizable=True, title='Match Pairs')    
    
    def _calculate_button_fired(self):
        self.calculate()

    def _save_button_fired(self):
        self.save_file()

    def _cluster_button_fired(self):
        self.cluster()

    def pick_handler(self, idx):
        txt = '%f %f %f    %f %f %f   %f %f %f' % tuple(self.pm.aPAIRS[idx][2].flat)
        self.txt = '%d:  %s' % (idx, txt) 
        
    def mk_v(self):
        spheres1 = gr.Spheres(phi_resolution=4, theta_resolution=4, input=self.p.poly_data, color_by='diffs', radius=0.001)
        spheres1.mapper.use_lookup_table_scalar_range = True
        spheres1.lut.range = (0, self.max_diff)
        spheres1.tv_glyph.generate_point_ids = True
        spheres1.tv_glyph.modified()
        v = ivtk.IVTK(size=(800,800))
        v.open()
        v.scene._renwin.stereo_type = 'anaglyph'
        v.scene._renwin.anaglyph_color_mask = 4,2
        v.scene.add_actor(spheres1.actor)
        self.v = v
        self.spheres = spheres1
        picker = Picker(v.scene, spheres1, self.pick_handler)
        v.scene.add_trait('picker', Instance(Picker))
        v.scene.picker = picker

    def calculate(self):
        self.UsAreAvailable = False 
        expt = self.ecf.h_unit_g
        mono = self.rf.rays.k_unit
        print '\nCalculating.... #expt=%d  #mono=%d\n' % (len(expt), len(mono))
	f= open("PATH_SIZE","w")
	f.write(str(len(expt)))
	f.close()
        #self.pm = laue_util.indexing.pairs_matching.PairsMatching(expt, mono, self.angular_cutoff)
        #self.pm.calculate()
        if self.show_3D:
            self.p = PointsAndBonds(points=self.pm.eulers)
            self.p.add_atom_data(1.0*self.pm.diffs, 'diffs')
            self.max_diff = np.max(self.pm.diffs)
            if self.v is None: self.mk_v()
            self.spheres.input = self.p.poly_data

    def cluster(self):
        self.c1.xyz = self.pm.eulers
        TTS = reduce(lambda x,y: x+y, self.c1.TTSs)
        fix = lambda a: np.average(a, axis=0)
        eulers = np.array([fix(np.take(self.pm.eulers, tt, axis=0)) for tt in TTS])
        Us = [transformations.euler_matrix(*eu)[:3,:3] for eu in eulers]
        lens = np.array([len(tt) for tt in TTS])
        ord = np.argsort(lens)[::-1]
        self.eulers = np.take(eulers, ord, axis=0)
        self.Us = np.take(Us, ord, axis=0)
        self.UsAreAvailable = True 
        self.lens = np.take(lens, ord)
        NF = lambda n: 0.5 * (1 + np.sqrt(1 + 8 * n))  # Definition
        # Arrays of (expt cluster - mono ray indices) for each cluster
        self.TTS = [TTS[i] for i in ord]
        print
        for i in xrange(min(len(self.Us), 30)):
            if (i == 0):
		file_object  = open("xxx.txt", "a")
		file_object.write('# size: '+str(self.lens[i])+' ==> '+str(NF(self.lens[i]))+',  eulers: '+str(self.eulers[i])+')')
                file_object.write('\n')
		file_object.close()

            print '# size: %d ==> %d,  eulers: %s)' % (self.lens[i], NF(self.lens[i]), self.eulers[i])
            print format_inp(self.crystal, self.Us[i])

    def save_file(self):
        if (self.UsAreAvailable == False):
           print 'Potential matrices must be computed first'
           return

        if (self.orient_matrix_selection is None):
           print 'One orientation matrix must be selected'
           return
        else:
           if not (self.orient_matrix_selection in xrange(self.Us.shape[0])):
              print 'Orientation matrix index is not in selection range '
              return

        if (self.infos_filename is None):
           print 'Update expt filename must be specified'
           return
        # Create a copy of the input expt hdf5 file
        infosfile = h5py.File(self.infos_filename, 'w')
        # Add information about the selected transformation matrix
        # obtained by average normalized H-vectors clustering
        infosfile_matrix = infosfile.create_group("Global_Rotation_Matrix")
        InitialMatrix = self.Us[self.orient_matrix_selection]

        # monochromatic info
        monocell = np.zeros((3,3))
        monocell[0,:] = np.array(fm['cell/A'])
        monocell[1,:] = np.array(fm['cell/B'])
        monocell[2,:] = np.array(fm['cell/C'])
        MonochromOrientMatrix = np.transpose(monocell)
        cell_parameters, Eulers, DirectOrientMatrix, DirectRotMatrix, ShapeMatrix = laue_util.indexing.find_transformation_mtx.information(InitialMatrix, MonochromOrientMatrix)

        # General information
        infosfile_matrix.create_dataset('Matrix_Status', data='Clustering')
        infosfile_matrix.create_dataset('Indirect_Rotation_Matrix', data=InitialMatrix)
        infosfile_matrix.create_dataset('Direct_Rotation_Matrix', data=DirectRotMatrix)
        infosfile_matrix.create_dataset('Direct_Orientation_Matrix', data=DirectOrientMatrix)
        infosfile_matrix.create_dataset('Cell_Parameters', data=cell_parameters)
        infosfile_matrix.create_dataset('Euler_Angles', data=Eulers)

        # Clustering information
        Clustering_Data_Group  = infosfile_matrix.create_group('Clustering_Data')
        Clustering_Data_Group.create_dataset('Euler_Angle_Cluster_Sizes', data=self.lens[self.orient_matrix_selection])

        # Refinement of global rotation matrix using (expts_clusters, mono_rays) pairs
        cluster = self.TTS[self.orient_matrix_selection]
        angle_pairs = [self.pm.PAIRS[i] for i in cluster]
        mono_rays_pair_indices = np.array([np.array(pair[1]) for pair in angle_pairs])
        expt_clusters_pair_indices = np.array([np.array(pair[0]) for pair in angle_pairs])

        mono_rays =  self.rf.rays.k_unit[np.ravel(mono_rays_pair_indices),:]
        expt_clusters = self.ecf.h_unit_g[np.ravel(expt_clusters_pair_indices),:]

        # Unique monorays indices and vectors
        Clustering_Data_Group.create_dataset('Mono_Rays', data=mono_rays)

        Clustering_Data_Group.create_dataset('Expt_Clusters', data=expt_clusters)

        # Infos files
        Input_Filenames_Group  = infosfile_matrix.create_group("Input_Filenames")
        Input_Filenames_Group.create_dataset('Used_Monochromatic_File', data='mono')
        Input_Filenames_Group.create_dataset('Used_Expt_File', data=sys.argv[2])

        infosfile.close()
        print 'Infofile written.'

a = fm['cell']['a'][...]
b = fm['cell']['b'][...]
c = fm['cell']['c'][...]
alpha = fm['cell']['alpha'][...]
beta = fm['cell']['beta'][...]
gamma = fm['cell']['gamma'][...]

m = Main()
m.ecf = ecf
m.rf = rf
m.edit_traits()
m.c1 = laue_util.cluster_points.Cluster1()
m.c1.paic_cell_size = 0.005
m.c1.link_cutoff = 0.005
m.c1.fcluster_cutoff = 0.005
m.c1.min_n = 10
c1_euler_ui = m.c1.edit_traits()
c1_euler_ui.title = 'Cluster eulers'          

m.crystal = (a, b, c, alpha, beta, gamma) 
m.group_no = 1

m.calculate()
#m.cluster()
#m.save_file()

print 'End'
import os
import time
time.sleep(1)
pid = os.getpid()
os.kill(pid, 9)





