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

import scipy.linalg
import numpy as np
import matplotlib.pyplot as plt

from enthought.traits.api import HasTraits, Instance, Int, Float, Button, Str, Bool, List, Array, \
                                 Range, Property, DelegatesTo, Trait
from enthought.traits.ui.api import View, Item, Group, TextEditor, RangeEditor

import transformations
import laue_util
import h5py, copy

import jdVis.gr as gr
# -------------------------------------------------------------------------------------------------
Nb_arguments = len(sys.argv)-1
# Nb_arguments different to 3 or 4
if ((Nb_arguments < 3) or (Nb_arguments > 4)):
    msg = '''
***************************************************************************************************
    
    Usage:
    
    ipython -wthread lu__refine_orientmatrix.py file_expt_extend.h5 file_mono.h5 file.mda (OPTIONAL infos.h5)
    
***************************************************************************************************
    '''
    sys.exit(msg)

# -------------------------------------------------------------------------------------------------
# Expt data
fe, de = laue_util.snop.read_snops(sys.argv[1])
expt = de['expt']
Totalnbframes = fe['lauecollect/frame_number_global'].shape[0]
FrameIndices = list(fe['lauecollect/frame_number_global'])
# Extract expt filename
fe_name = sys.argv[1].split('.')[0]

# Monochromatic data
fm, dm = laue_util.snop.read_snops(sys.argv[2])
rays = dm['rays_x']

# the mda file is required only to obtain the minimal lambda value
file_mda = sys.argv[3]
lambda_cutoff = 0.01
lambda_N = 120
mda = laue_util.read_mda(file_mda)                           # a dictionary of raw data from mda file
mda_f = laue_util.mda2lc.MdaToLC()
mda_f.lambda_cutoff = lambda_cutoff
mda_f.max_N = lambda_N
mda_f.mda_in = mda
lambda_min = min(mda_f.lc.l)

if (Nb_arguments == 4):
   infosfile_name = sys.argv[4]
# -------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------

class Main(HasTraits):

    # Graphical interface parameters

    # Expt Intensity Filter parameters
    intensity = Trait('I', 'I_bpf')
    # Cutoff specified by the user
    min_I = Float
    max_sintheta = Float

    # Initial Transformation Matrix selection
    manual_definition = Str
    manual_InitialGlobalMatrix = Array
    matrix_format = Trait('LaueUtil','LaueGUI')
    
    # Refinement parameters
    crystal_motion_max_angle = Float
    cell_shape_max_angle = Float
    min_nb_expected_pairs = Int
    max_nb_cycles = Int
   
    # Button 
    use_exptfile_matrix = Bool
    optimize_button = Button('Optimize')
    full_optimize_button = Button('Full optimize')
    use_exptfile_optimize_data = Bool
    file_name = Str
    infofile_button = Button('Save exptfile')

    # Comment window
    msg = Str

    # Internal variables
    opt_matrix = Instance(laue_util.optimize_TMatrix.Optimize_TMatrix)
    Totalnbframes = Int(0)
    FrameIndices = Array
    exptfile_InitialGlobalMatrix_evaluable = Bool
    exptfile_InitialGlobalMatrix = Array
    exptfile_PreviousRotationMatrices_evaluable = Bool
    exptfile_PreviousRotationMatrices = List

    # This variable is used to know the last optimization performance
    # It is used during output expt file writting
    LastOptimization = Trait('None','Manual2Rotation','Expt2Rotation','MRot2Trans','ERot2Trans','PRot2Trans')
    
    # Mono chromatic data
    rays = Instance(laue_util.Snop)
    # Laue data
    expt = Instance(laue_util.expt.Expt)


    # Selection of frame subset
    use_frame_selection = Bool
    selection_str = Str
    selection_list = List

    main_group = Group(
                       Group(
                               Item('intensity'),
                               Item('min_I', editor=TextEditor(evaluate=float, auto_set=False, enter_set=True)),
                               Item('max_sintheta', editor=TextEditor(evaluate=float, auto_set=False, enter_set=True))
                       ),
                       Group(
                               Item('use_exptfile_matrix',enabled_when='exptfile_InitialGlobalMatrix_evaluable'),
                               Item('manual_definition', enabled_when='not(use_exptfile_matrix)',editor=TextEditor(evaluate=str, auto_set=False, enter_set=True)),
                               Item('matrix_format',enabled_when='not(use_exptfile_matrix)')
                       ),
                       Group(
                               Item('crystal_motion_max_angle', editor=TextEditor(evaluate=float, auto_set=False, enter_set=True)),
                               Item('cell_shape_max_angle', editor=TextEditor(evaluate=float, auto_set=False, enter_set=True)),
                               Item('min_nb_expected_pairs', editor=TextEditor(evaluate=int, auto_set=False, enter_set=True)),
                               Item('max_nb_cycles_per_frame', editor=TextEditor(evaluate=int, auto_set=False, enter_set=True)),
                               Item('use_frame_selection'),
                               Item('selection_str', enabled_when='(use_frame_selection)',editor=TextEditor(evaluate=str, auto_set=False, enter_set=True))
                       ),
                       Group(
                               Item('optimize_button', show_label=False),
                               Item('full_optimize_button', show_label=False),
                               Item('use_exptfile_optimize_data',enabled_when='exptfile_PreviousRotationMatrices_evaluable')
                       ),
                       Group(
                               Item('file_name', editor=TextEditor(evaluate=str, auto_set=False, enter_set=True)),
                               Item('infofile_button', show_label=False)
                       ))
    traits_view = View(main_group, resizable=True, title='Orientation matrix optimization')    

    def _manual_definition_changed(self):
        self.manual_InitialGlobalMatrix = np.array(map(float, self.manual_definition.split())).reshape((3,3))

    def _infofile_button_fired(self):
        self.save_exptfile()

    def _optimize_button_fired(self):
        self.optimize()

    def _full_optimize_button_fired(self):
        self.full_optimize()

    def optimize(self):
        if (self.use_exptfile_matrix == True):
           InitialMatrix = np.copy(self.exptfile_InitialGlobalMatrix)
           self.opt_matrix.Used_InitialRU_Status = 'Clustering'
        else:
           InitialMatrix = np.copy(self.manual_InitialGlobalMatrix)
           self.opt_matrix.Used_InitialRU_Status = 'Manual'
           if (self.matrix_format == 'LaueGUI'):
              # conversion of coordinate systems to get a LU matrix from PG matrix
              Inv_CC = np.array([[0,0,1], [1,0,0], [0,1,0]])
              InitialMatrix = np.dot(Inv_CC, InitialMatrix)                    

        # Subset selection
        if self.use_frame_selection == True:
           self.process_selection()
        else:
            self.selection_list = FrameIndices[:]

        self.opt_matrix.filtering_expt(self.intensity, self.min_I, self.max_sintheta)
        # Optimization per frame of the rotation matrix
        self.opt_matrix.optimize(InitialMatrix, self.cell_shape_max_angle, self.crystal_motion_max_angle, self.max_nb_cycles_per_frame, self.min_nb_expected_pairs, self.selection_list, self.selection_str)

        # Summarize of the processing
        print "Number of frame rotation matrices not refined: %d / %d" % (self.opt_matrix.nbFailedRotMatConverge, self.opt_matrix.Totalnbframes)
        print "Number of frame rotation matrices not converged: %d / %d" % (self.opt_matrix.nbFailedRotMatRefine, self.opt_matrix.Totalnbframes)
        print self.FrameIndices[np.nonzero(self.opt_matrix.Array_FrameRU_converged == False)[0]]
        # Updating of the current optimization status
        if (self.use_exptfile_matrix == True):
           self.LastOptimization = 'Expt2Rotation'
        else:
           self.LastOptimization = 'Manual2Rotation'
        self.file_name = fe_name + '_rot.h5'

    def full_optimize(self):
        if (self.use_exptfile_optimize_data == True):
           self.opt_matrix.Used_List_initialRU_status = 'Previous'
           initial_ListRotationMatrix = self.exptfile_PreviousRotationMatrices
        else:
           if (self.opt_matrix.List_FrameRU is None): 
              print "Frame rotation matrices must be refined first"
              return
           self.opt_matrix.Used_List_initialRU_status = 'Refined'
           initial_ListRotationMatrix = self.opt_matrix.List_FrameRU 

        # Subset selection
        if self.use_frame_selection == True:
           self.process_selection()
        else:
           self.selection_list = FrameIndices[:]

        self.opt_matrix.filtering_expt(self.intensity, self.min_I, self.max_sintheta)
        # Optimization per frame of the transformation matrix
        self.opt_matrix.full_optimize(initial_ListRotationMatrix, self.cell_shape_max_angle, self.crystal_motion_max_angle, self.max_nb_cycles_per_frame, self.min_nb_expected_pairs, self.selection_list, self.selection_str)

        # Summarize of the processing
        print "Number of frame transformation matrices not refined: %d / %d" % (self.opt_matrix.nbFailedTransMatConverge, self.opt_matrix.Totalnbframes)
        print "Number of frame transformation matrices not converged: %d / %d" % (self.opt_matrix.nbFailedTransMatRefine, self.opt_matrix.Totalnbframes)
        print self.FrameIndices[np.nonzero(self.opt_matrix.Array_FrameTU_converged == False)[0]]

        # Updating of the current optimization status
        if (self.use_exptfile_optimize_data == False):
           if ((self.LastOptimization == 'Manual2Rotation') or (self.LastOptimization == 'MRot2Trans')): LastOptimization = 'MRot2Trans'
           if ((self.LastOptimization == 'Expt2Rotation') or (self.LastOptimization == 'ERot2Trans')): LastOptimization = 'ERot2Trans'
        else:
           LastOptimization = 'PRot2Trans' # Previous matrix
        self.LastOptimization = LastOptimization
        self.file_name = fe_name + '_full.h5'

    def save_exptfile(self):
        if (self.file_name is None):
           print 'Info filename must be specified'
           return

        if (self.LastOptimization == 'None'):
           print 'Nothing to save'
           return

        # Create an extra information hdf5 file
        new_infosfile = h5py.File(self.file_name, 'w')

        # Monochromatic info
        monocell = np.zeros((3,3))
        monocell[0,:] = np.array(fm['cell/A'])
        monocell[1,:] = np.array(fm['cell/B'])
        monocell[2,:] = np.array(fm['cell/C'])
        MonochromOrientMatrix = np.transpose(monocell)

        # Global rotation matrix
        if ((self.LastOptimization == 'Expt2Rotation') or (self.LastOptimization == 'ERot2Trans') or (self.LastOptimization == 'PRot2Trans')):
           infosfile = h5py.File(infosfile_name, 'r')
           #new_infosfile_globalrotmatrix = new_infosfile.create_group("Global_Rotation_Matrix")

           # Expt Global Rotation Matrix
           infosfile.copy('Global_Rotation_Matrix',new_infosfile)
           
           # Previous Rotation Matrix
           if (self.LastOptimization == 'PRot2Trans'):
              #new_infosfile_framerotmatrices = new_infosfile.create_group("Frame_Rotation_Matrices")
              infosfile.copy('Frame_Rotation_Matrices',new_infosfile)

           infosfile.close()

        # New data
        if ((self.LastOptimization == 'Manual2Rotation') or (self.LastOptimization == 'MRot2Trans')):
           # Global Rotation Matrix
           new_infosfile_globalrotmatrix = new_infosfile.create_group("Global_Rotation_Matrix")
           cell_parameters, Eulers, DirectOrientMatrix, DirectRotMatrix, ShapeMatrix = laue_util.indexing.find_transformation_mtx.information(self.opt_matrix.Used_InitialRU, MonochromOrientMatrix)
           InitialMatrix = self.opt_matrix.Used_InitialRU
           InitialMatrix_status = 'Manual'
           new_infosfile_globalrotmatrix.create_dataset('Indirect_Orientation_Matrix', data=self.opt_matrix.Used_InitialRU)
           new_infosfile_globalrotmatrix.create_dataset('Direct_Orientation_Matrix', data=DirectOrientMatrix)
           new_infosfile_globalrotmatrix.create_dataset('Direct_Rotation_Matrix', data=DirectRotMatrix)
           new_infosfile_globalrotmatrix.create_dataset('Matrix_Status', data=InitialMatrix_status)
           new_infosfile_globalrotmatrix.create_dataset('Cell_Parameters', data=cell_parameters)
           new_infosfile_globalrotmatrix.create_dataset('Euler_Angles', data=Eulers)
           new_infosfile_globalrotmatrix.create_dataset('Indirect_Rotation_Matrix', data=np.transpose(np.linalg.inv(DirectRotMatrix))) # correction

        # Frame Rotation matrix
        if ((self.LastOptimization == 'Manual2Rotation') or (self.LastOptimization == 'Expt2Rotation') or (self.LastOptimization == 'MRot2Trans') or (self.LastOptimization == 'ERot2Trans')):
           new_infosfile_matrix = new_infosfile.create_group("Frame_Rotation_Matrices")
           Indirect_Rotation_Matrices = []
           Direct_Rotation_Matrices = []
           Direct_Orientation_Matrices = []
           Cell_Parameters = []
           Euler_Angles = []

           for matrix in self.opt_matrix.List_FrameRU:
              Indirect_Rotation_Matrices.append(np.ravel(matrix, order='C'))
              cell_parameters, Eulers, DirectOrientMatrix, DirectRotMatrix, ShapeMatrix = laue_util.indexing.find_transformation_mtx.information(matrix, MonochromOrientMatrix)
              Direct_Rotation_Matrices.append(np.ravel(DirectRotMatrix, order='C'))
              Direct_Orientation_Matrices.append(np.ravel(DirectOrientMatrix, order='C'))
              Cell_Parameters.append(cell_parameters)
              Euler_Angles.append(Eulers)

           FrameRotMatrices_succeed = self.opt_matrix.Array_FrameRU_succeed
           new_infosfile_matrix.create_dataset('Indirect_Rotation_Matrices', data=Indirect_Rotation_Matrices)
           new_infosfile_matrix.create_dataset('Direct_Rotation_Matrices', data=Direct_Rotation_Matrices)
           new_infosfile_matrix.create_dataset('Refinement_Succeed', data=FrameRotMatrices_succeed)
           new_infosfile_matrix.create_dataset('Direct_Orientation_Matrices', data=Direct_Orientation_Matrices)
           new_infosfile_matrix.create_dataset('Cell_Parameters', data=Cell_Parameters)
           new_infosfile_matrix.create_dataset('Euler_Angles', data=Euler_Angles)
           new_infosfile_refine_parameters = new_infosfile_matrix.create_group("Refinement_Parameters")
           new_infosfile_refine_parameters.create_dataset('Crystal_Motion_Angle', data=self.opt_matrix.Array_FrameR_crystal_motion_max_angle)
           new_infosfile_refine_parameters.create_dataset('Crystal_Shape_Angle', data=self.opt_matrix.Array_FrameR_crystal_shape_max_angle)
           new_infosfile_refine_parameters.create_dataset('Min_I', data=self.opt_matrix.Array_FrameR_min_I)
           new_infosfile_refine_parameters.create_dataset('Nb_Pairs', data=self.opt_matrix.Array_FrameR_nbpairs)
           new_infosfile_initial = new_infosfile_refine_parameters.create_group("Initial_Parameters")
           for key in self.opt_matrix.R_arguments.keys():
              new_infosfile_initial.attrs[key] = self.opt_matrix.R_arguments[key]
           Array_FrameR_indices = np.array([], dtype=int)
           for i, FrameR_indices in enumerate(self.opt_matrix.List_FrameR_Indices):
              # Test if the current frame has pairs of indices
              if (FrameR_indices.shape[0] == 0): continue
              if (Array_FrameR_indices.shape[0] == 0):
                 Array_FrameR_indices = np.concatenate((i*np.ones(FrameR_indices.shape[0], dtype=int)[:,np.newaxis], FrameR_indices),axis=1)
              else:
                 Array_FrameR_indices = np.concatenate((Array_FrameR_indices, np.concatenate((i*np.ones(FrameR_indices.shape[0], dtype=int)[:,np.newaxis], FrameR_indices),axis=1)), axis=0)
           new_infosfile_refine_parameters.create_dataset('Pair_indices', data=Array_FrameR_indices)

        # Frame Transformation matrix
        if ((self.LastOptimization == 'MRot2Trans') or (self.LastOptimization == 'ERot2Trans') or (self.LastOptimization == 'PRot2Trans')):
             
           new_infosfile_matrix = new_infosfile.create_group("Frame_Transformation_Matrices")
           Indirect_Transformation_Matrices = []
           Direct_Rotation_Matrices = []
           Direct_Shape_Matrices = []
           Direct_Orientation_Matrices = []
           Cell_Parameters = []
           Euler_Angles = []

           for matrix in self.opt_matrix.List_FrameTU:
              Indirect_Transformation_Matrices.append(np.ravel(matrix, order='C'))
              cell_parameters, Eulers, DirectOrientMatrix, DirectRotMatrix, ShapeMatrix = laue_util.indexing.find_transformation_mtx.information(matrix, MonochromOrientMatrix)
              Direct_Rotation_Matrices.append(np.ravel(DirectRotMatrix, order='C'))
              Direct_Shape_Matrices.append(np.ravel(ShapeMatrix, order='C'))
              Direct_Orientation_Matrices.append(np.ravel(DirectOrientMatrix, order='C'))
              Cell_Parameters.append(cell_parameters)
              Euler_Angles.append(Eulers)

           FrameTransMatrices_succeed = self.opt_matrix.Array_FrameTU_succeed
           new_infosfile_matrix.create_dataset('Indirect_Transformation_Matrices', data=Indirect_Transformation_Matrices)
           new_infosfile_matrix.create_dataset('Direct_Rotation_Matrices', data=Direct_Rotation_Matrices)
           new_infosfile_matrix.create_dataset('Refinement_Succeed', data=FrameTransMatrices_succeed)
           new_infosfile_matrix.create_dataset('Direct_Orientation_Matrices', data=Direct_Orientation_Matrices)
           new_infosfile_matrix.create_dataset('Direct_Shape_Matrices', data=Direct_Shape_Matrices)
           new_infosfile_matrix.create_dataset('Cell_Parameters', data=Cell_Parameters)
           new_infosfile_matrix.create_dataset('Euler_Angles', data=Euler_Angles)
           new_infosfile_refine_parameters = new_infosfile_matrix.create_group("Refinement_Parameters")
           new_infosfile_refine_parameters.create_dataset('Crystal_Motion_Angle', data=self.opt_matrix.Array_FrameT_crystal_motion_max_angle)
           new_infosfile_refine_parameters.create_dataset('Crystal_Shape_Angle', data=self.opt_matrix.Array_FrameT_crystal_shape_max_angle)
           new_infosfile_refine_parameters.create_dataset('Min_I', data=self.opt_matrix.Array_FrameT_min_I)
           new_infosfile_refine_parameters.create_dataset('Nb_Pairs', data=self.opt_matrix.Array_FrameT_nbpairs)
           new_infosfile_initial = new_infosfile_refine_parameters.create_group("Initial_Parameters")
           for key in self.opt_matrix.T_arguments.keys():
              new_infosfile_initial.attrs[key] = self.opt_matrix.T_arguments[key]
           Array_FrameT_indices = np.array([], dtype=int)
           for i, FrameT_indices in enumerate(self.opt_matrix.List_FrameT_Indices):
              if (FrameT_indices.shape[0] == 0): continue
              if (Array_FrameT_indices.shape[0] == 0):
                 Array_FrameT_indices = np.concatenate((i*np.ones(FrameT_indices.shape[0], dtype=int)[:,np.newaxis], FrameT_indices),axis=1)
              else:
                 Array_FrameT_indices = np.concatenate((Array_FrameT_indices, np.concatenate((i*np.ones(FrameT_indices.shape[0], dtype=int)[:,np.newaxis], FrameT_indices),axis=1)), axis=0)
           new_infosfile_refine_parameters.create_dataset('Pair_indices', data=Array_FrameT_indices)
        new_infosfile.close()
        print 'New infosfile written.'

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
           
m = Main()
# Initialize parameter values
m.intensity = 'I'
m.min_I = 1000
m.max_sintheta = 0.0
# Initial Transformation Matrix selection
m.use_exptfile_matrix = True
m.manual_definition = ''
m.matrix_format = 'LaueGUI'
# Refinement parameters
m.crystal_motion_max_angle = 2.0
m.cell_shape_max_angle = 1.0
m.min_nb_expected_pairs = 20
m.max_nb_cycles_per_frame = 5
# Button 
m.file_name = ''
m.use_frame_selection = False
m.selection_str = ''

if (Nb_arguments == 4):
   infosfile = h5py.File(sys.argv[4], 'r')

   # Check if initial matrix evaluable in expt file
   if ('Global_Rotation_Matrix' in infosfile.listnames()):
      m.exptfile_InitialGlobalMatrix_evaluable = True
      m.exptfile_InitialGlobalMatrix = np.array(infosfile['Global_Rotation_Matrix/Indirect_Rotation_Matrix'])
      m.use_exptfile_matrix = True
   else:
      m.exptfile_InitialGlobalMatrix_evaluable = False
      m.use_exptfile_matrix = False

   # Check if frame rotation matrices evaluable in expt file
   if ('Frame_Rotation_Matrices' in infosfile.listnames()):
      m.exptfile_PreviousRotationMatrices_evaluable = True
      nb_frames = infosfile['Frame_Rotation_Matrices/Indirect_Rotation_Matrices'].shape[0]
      m.exptfile_PreviousRotationMatrices = [np.array(infosfile['Frame_Rotation_Matrices/Indirect_Rotation_Matrices'][k,:]).reshape((3,3)) for k in xrange(nb_frames)]
      m.use_exptfile_optimize_data = True
   else:
      m.exptfile_PreviousRotationMatrices_evaluable = False
      m.use_exptfile_optimize_data = False
   infosfile.close()
else:
      m.exptfile_InitialGlobalMatrix_evaluable = False
      m.use_exptfile_matrix = False
      m.exptfile_PreviousRotationMatrices_evaluable = False
      m.use_exptfile_optimize_data = False
 
m.rays = rays
m.expt = expt
m.LastOptimization = 'None'
m.Totalnbframes = Totalnbframes
m.FrameIndices = np.array(FrameIndices)
# orientation matrix optimization object
m.opt_matrix = laue_util.optimize_TMatrix.Optimize_TMatrix(expt,rays,lambda_min,Totalnbframes,FrameIndices)
m.edit_traits()

m.optimize()
m.save_exptfile()

print 'End'
import os
import time
time.sleep(1)
pid = os.getpid()
os.kill(pid, 9)
 

