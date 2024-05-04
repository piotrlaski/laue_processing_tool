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
import sys

import numpy as np
from scipy import stats, optimize
import h5py

import wx
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure
from enthought.traits.api import HasTraits, Instance, Button, Str, Int, Float, Bool, Trait, Array
from enthought.traits.ui.api import View, Item, Group, CustomEditor, TextEditor, ColorEditor

import laue_util
from laue_util.ui.sfp import _mk_fig, SingleFigurePanel

import string
import transformations


# -------------------------------------------------------------------------------------------------

if not(len(sys.argv) == 2):
    msg = '''
***************************************************************************************************
    
    Usage:
    
    ipython -wthread lu__view_orientmatrix.py infos.h5
    
***************************************************************************************************
    '''
    sys.exit(msg)

# -------------------------------------------------------------------------------------------------

class MAIN(SingleFigurePanel):

    nb_frames = Int
    Plot_menu = Trait('Euler_angles','Cell_angles','Cell_lengths','Nb_pairs')
    Global_rotation = Bool
    Frame_rotation = Bool
    Frame_transformation = Bool

    # Local variables
    Global_Rotation_Matrix_evaluable = Bool
    Global_Rotation_Cell_Parameters = Array
    Global_Rotation_Euler_Angles = Array
    Global_Rotation_Matrix = Array

    Frame_Rotation_Matrix_evaluable = Bool
    Frame_Rotation_Cell_Parameters = Array
    Frame_Rotation_Euler_Angles = Array
    Frame_Rotation_Matrices = Array
    Shift_Frame_Rotation_Eulers = Array
    Max_Frame_Rotation_Angle = Array

    Frame_Rotation_Nb_Used_Pairs = Array

    Frame_Transformation_Matrix_evaluable = Bool
    Frame_Transformation_Cell_Parameters = Array
    Frame_Transformation_Euler_Angles = Array
    Frame_Transformation_Matrices = Array
    Shift_Frame_Transformation_Eulers = Array
    Max_Frame_Transformation_Angle = Array

    Frame_Transformation_Nb_Used_Pairs = Array

    Frame_indices = Array

    Global_Rotation_Matrix_option_evaluable = Bool
    Frame_Rotation_Matrix_option_evaluable = Bool
    Frame_Transformation_Matrix_option_evaluable = Bool

    Frame_Rotation_Nb_Pairs_evaluable = Bool
    Frame_Transformation_Nb_Pairs_evaluable = Bool

    legend = Str
    colors = Str
    lines = Str
    # Graphical interface
    traits_view = View(Group(
                             Group( 
                                 Item('Plot_menu'),
                                 Item('Global_rotation', enabled_when='Global_Rotation_Matrix_option_evaluable'),
                                 Item('Frame_rotation', enabled_when='Frame_Rotation_Matrix_option_evaluable'),
                                 Item('Frame_transformation', enabled_when='Frame_Transformation_Matrix_option_evaluable'),
                                 Item(name='legend', style='readonly', editor=TextEditor(multi_line=True))
                                 ),
                             Group(
                                   Item('fig_mono', editor=CustomEditor(_mk_fig),
                                        resizable=True, 
                                        show_label=False),
                             ),
                             orientation='vertical'
                             ),
                       resizable=True,
                       title = 'Matrix_Parameters_Viewer'
                       )
    
    def __init__(self, f1, *args, **kw):

        infosfile = h5py.File(f1, 'r')

        # Check if initial matrix evaluable
        self.Global_Rotation_Matrix_evaluable = True
        self.Global_Rotation_Cell_Parameters = np.array(infosfile['Global_Rotation_Matrix/Cell_Parameters'])
        self.Global_Rotation_Euler_Angles = np.array(infosfile['Global_Rotation_Matrix/Euler_Angles'])
        self.Global_Rotation_Matrix = np.array(infosfile['Global_Rotation_Matrix/Direct_Rotation_Matrix'])
        Inv_Global_Rotation_Matrix = np.linalg.inv(self.Global_Rotation_Matrix)
        # Check if frame rotation matrices evaluable
        if ('Frame_Rotation_Matrices' in infosfile.listnames()):
           self.Frame_Rotation_Matrix_evaluable = True
           self.Frame_Rotation_Cell_Parameters = np.transpose(infosfile['Frame_Rotation_Matrices/Cell_Parameters'])
           self.Frame_Rotation_Euler_Angles = np.transpose(infosfile['Frame_Rotation_Matrices/Euler_Angles'])
           self.Nb_Frames = self.Frame_Rotation_Euler_Angles.shape[1]
           self.Frame_indices = np.arange(self.Nb_Frames)
           Frame_Rotation_Matrices = np.array(infosfile['Frame_Rotation_Matrices/Direct_Rotation_Matrices'])
           self.Frame_Rotation_Matrices = [Frame_Rotation_Matrices[k,:].reshape((3,3)) for k in xrange(self.Nb_Frames)]
           Shift_Frame_Rotation_Matrices = [np.dot(self.Frame_Rotation_Matrices[k],Inv_Global_Rotation_Matrix) for k in xrange(self.Nb_Frames)]
           Shift_Frame_Rotation_Eulers = [transformations.euler_from_matrix(Shift_Frame_Rotation_Matrices[k]) for k in xrange(self.Nb_Frames)]
           self.Shift_Frame_Rotation_Eulers = np.rad2deg(np.transpose(np.array(Shift_Frame_Rotation_Eulers)))

           if ('Nb_Pairs' in infosfile['Frame_Rotation_Matrices/Refinement_Parameters'].listnames()):
              self.Frame_Rotation_Nb_Used_Pairs = np.array(infosfile['Frame_Rotation_Matrices/Refinement_Parameters/Nb_Pairs'])
              self.Frame_Rotation_Nb_Pairs_evaluable = True
           else:
              self.Frame_Rotation_Nb_Pairs_evaluable = False
        else:
           self.Frame_Rotation_Matrix_evaluable = False

        # Check if frame rotation matrices evaluable
        if ('Frame_Transformation_Matrices' in infosfile.listnames()):
           self.Frame_Transformation_Matrix_evaluable = True
           self.Frame_Transformation_Cell_Parameters = np.transpose(infosfile['Frame_Transformation_Matrices/Cell_Parameters'])
           self.Frame_Transformation_Euler_Angles = np.transpose(infosfile['Frame_Transformation_Matrices/Euler_Angles'])

           Frame_Transformation_Matrices = np.array(infosfile['Frame_Transformation_Matrices/Direct_Rotation_Matrices'])
           self.Frame_Transformation_Matrices = [Frame_Transformation_Matrices[k,:].reshape((3,3)) for k in xrange(self.Nb_Frames)]
           Shift_Frame_Transformation_Matrices = [np.dot(self.Frame_Transformation_Matrices[k],Inv_Global_Rotation_Matrix) for k in xrange(self.Nb_Frames)]
           Shift_Frame_Transformation_Eulers = [transformations.euler_from_matrix(Shift_Frame_Transformation_Matrices[k]) for k in xrange(self.Nb_Frames)]
           self.Shift_Frame_Transformation_Eulers = np.rad2deg(np.transpose(np.array(Shift_Frame_Transformation_Eulers)))

           if ('Nb_Pairs' in infosfile['Frame_Transformation_Matrices/Refinement_Parameters'].listnames()):
              self.Frame_Transformation_Nb_Used_Pairs = np.array(infosfile['Frame_Transformation_Matrices/Refinement_Parameters/Nb_Pairs'])
              self.Frame_Transformation_Nb_Pairs_evaluable = True
           else:
              self.Frame_Transformation_Nb_Pairs_evaluable = False
        else:
           self.Frame_Transformation_Matrix_evaluable = False

        infosfile.close()
        if not((self.Frame_Rotation_Matrix_evaluable == True) or (self.Frame_Transformation_Matrix_evaluable == True)):
           print 'Only global rotation matrix available'
           return

        super(MAIN, self).__init__(*args, **kw)

    def _Global_rotation_changed(self):
        self.calculate()

    def _Frame_rotation_changed(self):
        self.calculate()

    def _Frame_transformation_changed(self):
        self.calculate()

    def _Plot_menu_changed(self):
        if (self.Plot_menu == 'Euler_angles'):
           self.Global_rotation = False
           self.Global_Rotation_Matrix_option_evaluable = False
           self.Frame_Rotation_Matrix_option_evaluable = self.Frame_Rotation_Matrix_evaluable
           self.Frame_Transformation_Matrix_option_evaluable = self.Frame_Transformation_Matrix_evaluable

        if (self.Plot_menu == 'Cell_angles'):
           self.Frame_rotation = False
           self.Global_Rotation_Matrix_option_evaluable = self.Global_Rotation_Matrix_evaluable
           self.Frame_Rotation_Matrix_option_evaluable = False
           self.Frame_Transformation_Matrix_option_evaluable = self.Frame_Transformation_Matrix_evaluable

        if (self.Plot_menu == 'Cell_lengths'):
           self.Frame_rotation = False
           self.Global_Rotation_Matrix_option_evaluable = self.Global_Rotation_Matrix_evaluable
           self.Frame_Rotation_Matrix_option_evaluable = False
           self.Frame_Transformation_Matrix_option_evaluable = self.Frame_Transformation_Matrix_evaluable

        if (self.Plot_menu == 'Nb_pairs'):
           self.Frame_rotation = False
           self.Global_Rotation_Matrix_option_evaluable = False
           self.Frame_Rotation_Matrix_option_evaluable = self.Frame_Rotation_Nb_Pairs_evaluable
           self.Frame_Transformation_Matrix_option_evaluable = self.Frame_Transformation_Nb_Pairs_evaluable

        self.calculate()

    def calculate(self):
        self.graph()
        self.calc()

    def graph(self):
        self.AxisX_label = str('Frame indices')

        if (self.Plot_menu == 'Euler_angles'):
           self.AxisY_label = str('Euler angles')
	   self.Title = str('Plot Euler angles - Frame indices')
           self.colors = ''
           self.colors += 'euler1 in red \n'
           self.colors += 'euler2 in green \n'
           self.colors += 'euler3 in blue \n'

        if (self.Plot_menu == 'Cell_angles'):
           self.AxisY_label = str('Cell angles')
	   self.Title = str('Plot Cell angles - Frame indices')
           self.colors = ''
           self.colors += 'alpha in red \n'
           self.colors += 'beta in green \n'
           self.colors += 'gamma in blue \n'

        if (self.Plot_menu == 'Cell_lengths'):
           self.AxisY_label = str('Cell lengths')
	   self.Title = str('Plot Cell lengths - Frame indices')
           self.colors = ''
           self.colors += 'a in red \n'
           self.colors += 'b in green \n'
           self.colors += 'c in blue \n'

        if (self.Plot_menu == 'Nb_pairs'):
           self.AxisY_label = str('Nb pairs')
	   self.Title = str('Plot Nb pairs - Frame indices')
           self.colors = ''

    def calc(self):

        self.clear()
        # Legend for lines
        self.lines = ''
        self.legend = ''
        # Plot of euler angles 
        if (self.Plot_menu == 'Euler_angles'):
           if (self.Frame_rotation == True):
              self.lines += 'Frame rotations in straight lines \n'

              print self.Nb_Frames
              self.ax.axis([0, self.Nb_Frames, -10, 10])
#              self.ax.xlim(0, self.Nb_Frames)
              self.plot = self.ax.plot(self.Frame_indices, self.Shift_Frame_Rotation_Eulers[0,:], linestyle='-', color='red', markersize=5)
              self.plot = self.ax.plot(self.Frame_indices, self.Shift_Frame_Rotation_Eulers[1,:], linestyle='-', color='green', markersize=5)
              self.plot = self.ax.plot(self.Frame_indices, self.Shift_Frame_Rotation_Eulers[2,:], linestyle='-', color='blue', markersize=5)

           if (self.Frame_transformation == True):
              self.lines += 'Frame transformations in dash lines \n'

              print self.Nb_Frames
              self.ax.axis([0, self.Nb_Frames, -10, 10])
#              self.ax.xlim(0, self.Nb_Frames)
              self.plot = self.ax.plot(self.Frame_indices, self.Shift_Frame_Transformation_Eulers[0,:], linestyle='--', color='red', markersize=5)
              self.plot = self.ax.plot(self.Frame_indices, self.Shift_Frame_Transformation_Eulers[1,:], linestyle='--', color='green', markersize=5)
              self.plot = self.ax.plot(self.Frame_indices, self.Shift_Frame_Transformation_Eulers[2,:], linestyle='--', color='blue', markersize=5)

        # Plot of cell angles 
        if (self.Plot_menu == 'Cell_angles'):
           if (self.Global_rotation == True):
              self.lines += 'Global rotation in straight lines \n'

              Global_Rotation = np.zeros((3,self.Nb_Frames))
              Global_Rotation[0,:] = self.Global_Rotation_Cell_Parameters[3]
              Global_Rotation[1,:] = self.Global_Rotation_Cell_Parameters[4]
              Global_Rotation[2,:] = self.Global_Rotation_Cell_Parameters[5]

              self.plot = self.ax.plot(self.Frame_indices, Global_Rotation[0,:], linestyle='-', color='red', markersize=5)
              self.plot = self.ax.plot(self.Frame_indices, Global_Rotation[1,:], linestyle='-', color='green', markersize=5)
              self.plot = self.ax.plot(self.Frame_indices, Global_Rotation[2,:], linestyle='-', color='blue', markersize=5)

           if (self.Frame_transformation == True):
              self.lines += 'Frame transformations in dash lines \n'

              self.plot = self.ax.plot(self.Frame_indices, self.Frame_Transformation_Cell_Parameters[3,:], linestyle='--', color='red', markersize=5)
              self.plot = self.ax.plot(self.Frame_indices, self.Frame_Transformation_Cell_Parameters[4,:], linestyle='--', color='green', markersize=5)
              self.plot = self.ax.plot(self.Frame_indices, self.Frame_Transformation_Cell_Parameters[5,:], linestyle='--', color='blue', markersize=5)

        # Plot of transformation angles 
        if (self.Plot_menu == 'Cell_lengths'):
           if (self.Global_rotation == True):
              self.lines += 'Global rotation in straight lines \n'

              Global_Rotation = np.zeros((3,self.Nb_Frames))
              Global_Rotation[0,:] = self.Global_Rotation_Cell_Parameters[0]
              Global_Rotation[1,:] = self.Global_Rotation_Cell_Parameters[1]
              Global_Rotation[2,:] = self.Global_Rotation_Cell_Parameters[2]

              self.plot = self.ax.plot(self.Frame_indices, Global_Rotation[0,:], linestyle='-', color='red', markersize=5)
              self.plot = self.ax.plot(self.Frame_indices, Global_Rotation[1,:], linestyle='-', color='green', markersize=5)
              self.plot = self.ax.plot(self.Frame_indices, Global_Rotation[2,:], linestyle='-', color='blue', markersize=5)

           if (self.Frame_transformation == True):
              self.lines += 'Frame transformations in dash lines \n'
              self.plot = self.ax.plot(self.Frame_indices, self.Frame_Transformation_Cell_Parameters[0,:], linestyle='--', color='red', markersize=5)
              self.plot = self.ax.plot(self.Frame_indices, self.Frame_Transformation_Cell_Parameters[1,:], linestyle='--', color='green', markersize=5)
              self.plot = self.ax.plot(self.Frame_indices, self.Frame_Transformation_Cell_Parameters[2,:], linestyle='--', color='blue', markersize=5)

        # Plot of nb pairs 
        if (self.Plot_menu == 'Nb_pairs'):
           if (self.Frame_rotation == True):
              self.lines += 'Frame rotations in straight lines \n'
              self.plot = self.ax.plot(self.Frame_indices, self.Frame_Rotation_Nb_Used_Pairs, linestyle='-', color='red', markersize=5)

           if (self.Frame_transformation == True):
              self.lines += 'Frame transformations in dash lines \n'
              self.plot = self.ax.plot(self.Frame_indices, self.Frame_Transformation_Nb_Used_Pairs, linestyle='--', color='red', markersize=5)

        if not(self.lines==''):
           self.legend = self.lines + '\n' + self.colors

        self.plot = self.ax.set_xlabel(self.AxisX_label)
        self.plot = self.ax.set_ylabel(self.AxisY_label)
        self.plot = self.ax.set_title(self.Title)
        self.ax.grid()
        self.draw()

        self.ax.figure.savefig('euler_plot.png')

# -------------------------------------------------------------------------------------------------

m = MAIN(sys.argv[1])
m.edit_traits()
m.Plot_menu = 'Euler_angles'
m.Global_Rotation_Matrix_option_evaluable = False
m.Frame_Rotation_Matrix_option_evaluable = True
m.Frame_Transformation_Matrix_option_evaluable = True
m.Global_rotation = False
m.Frame_rotation = True
m.Frame_transformation = False
m.calculate()

print 'End'
import os
import time
time.sleep(1)
pid = os.getpid()
os.kill(pid, 9)
 

