#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Feb 15, 2010 by: rch

from etsproxy.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Bool, Enum, Event, implements, DelegatesTo, \
    Callable

import os

from numpy import \
    argmax, frompyfunc

from mathkit.mfn import MFnLineArray

#-- Tabular Adapter Definition -------------------------------------------------
from os.path import exists

from etsproxy.traits.ui.api \
    import View, Item, VGroup, HGroup

from matresdev.db.exdb.ex_type import ExType
from matresdev.db.exdb.i_ex_type import IExType

from matresdev.db.matdb.trc.concrete_mixture \
    import ConcreteMixture
    
from matresdev.db.simdb import \
    SimDB

# Access to the toplevel directory of the database
#
simdb = SimDB()

#class ExpBendingTestThreePoint(ExType):
class ExpGrouting(ExType):
    '''Experiment: Bending Test Three Point 
    '''
#    label = Str('three point bending test')

    implements(IExType)

    file_ext = 'raw'

    #--------------------------------------------------------------------
    # register a change of the traits with metadata 'input'
    #--------------------------------------------------------------------

    input_change = Event
    @on_trait_change('+input, ccs.input_change, +ironing_param')
    def _set_input_change(self):
        self.input_change = True

    #--------------------------------------------------------------------------------
    # specify inputs:
    #--------------------------------------------------------------------------------

    length = Float(1.15, unit = 'm', input = True, table_field = True,
                           auto_set = False, enter_set = True)
    width = Float(0.2, unit = 'm', input = True, table_field = True,
                           auto_set = False, enter_set = True)
    thickness = Float(0.06, unit = 'm', input = True, table_field = True,
                           auto_set = False, enter_set = True)

    #--------------------------------------------------------------------------
    # composite cross section
    #--------------------------------------------------------------------------


    cm = Instance(ConcreteMixture)
    def _cm_default(self):
        '''default settings correspond to 
        setup 'PZ-0708-1'
        '''
        return ConcreteMixture.db['PZ-0708']

    #--------------------------------------------------------------------------------
    # define processing
    #--------------------------------------------------------------------------------

    def _read_data_array(self):
        ''' Read the experiment data. 
        '''
        if exists(self.data_file):

            print 'READ FILE'
            # change the file name dat with asc  
            file_split = self.data_file.split('.')

            file_name = file_split[0] + '.csv'
            if not os.path.exists(file_name):

                file_name = file_split[0] + '.raw'
                if not os.path.exists(file_name):
                    raise IOException, 'file %s does not exist' % file_name

            print 'file_name', file_name

            _data_array = loadtxt_bending(file_name)

            self.data_array = _data_array

    names_and_units = Property
    @cached_property
    def _get_names_and_units(self):
        ''' Set the names and units of the measured data.
        '''
        names = ['w', 'eps_c', 'F']
        units = ['mm', '%', 'kN']
        return names, units

#    mfn_elastomer = Instance( MFnLineArray )
#    def _mfn_elastomer_default( self ):
#        elastomer_path = os.path.join( simdb.exdata_dir, 'bending_tests', 'ZiE_2011-06-08_BT-12c-6cm-0-TU', 'elastomer_f-w.raw' )
#        # loadtxt_bending returns an array with three columns: 
#        # 0: deformation; 1: eps_c; 2: force
#        _data_array_elastomer = loadtxt_bending( elastomer_path )
#        return MFnLineArray( xdata = _data_array_elastomer[:, 0], ydata = _data_array_elastomer[:, 2] )
#
#    elastomer_force = Callable
#    def elastomer_force_default( self ):
#        return frompyfunc( self.mfn_elastomer.get_value, 2, 1 )

    def process_source_data(self):
        '''read in the measured data from file and assign
        attributes after array processing.        
        '''
        super(ExpBT3PT, self).process_source_data()


        elastomer_path = os.path.join(simdb.exdata_dir, 'bending_tests', 'three_point', '2011-06-10_BT-3PT-12c-6cm-0-TU_ZiE', 'elastomer_f-w.raw')
        _data_array_elastomer = loadtxt_bending(elastomer_path)

        # force [kN]:
        #
        xdata = -0.001 * _data_array_elastomer[:, 2].flatten()

        # displacement [mm]:
        #
        ydata = -1.0 * _data_array_elastomer[:, 0].flatten()

        mfn_displacement_elastomer = MFnLineArray(xdata = xdata, ydata = ydata)
        displacement_elastomer_vectorized = frompyfunc(mfn_displacement_elastomer.get_value, 1, 1)

        # convert data from 'N' to 'kN' and change sign
        #
        self.F = -0.001 * self.F

        # change sign in positive values for vertical displacement [mm]
        #
        self.w = -1.0 * self.w

        # substract the deformation of the elastomer cushion between the cylinder
        # 
        self.w = self.w - displacement_elastomer_vectorized(self.F)



    #--------------------------------------------------------------------------------
    # plot templates
    #--------------------------------------------------------------------------------

    plot_templates = {'force / deflection (center)'          : '_plot_force_deflection_center',
                      'smoothed force / deflection (center)' : '_plot_smoothed_force_deflection_center',
                     }

    default_plot_template = 'force / deflection (center)'

    def _plot_force_deflection_center(self, axes):
        xkey = 'deflection [mm]'
        ykey = 'force [kN]'
        # NOTE: processed data returns positive values for force and displacement
        #
        xdata = self.w
        ydata = self.F

        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))
        axes.plot(xdata, ydata
                       # color = c, linewidth = w, linestyle = s 
                       )

    n_fit_window_fraction = Float(0.1)

    def _plot_smoothed_force_deflection_center(self, axes):

        # get the index of the maximum stress
        max_force_idx = argmax(self.F)
        # get only the ascending branch of the response curve
        f_asc = self.F[:max_force_idx + 1]
        w_asc = self.w[:max_force_idx + 1]

        f_max = f_asc[-1]
        w_max = w_asc[-1]

        n_points = int(self.n_fit_window_fraction * len(w_asc))
        f_smooth = smooth(f_asc, n_points, 'flat')
        w_smooth = smooth(w_asc, n_points, 'flat')

        axes.plot(w_smooth, f_smooth, color = 'blue', linewidth = 2)

#        secant_stiffness_w10 = ( f_smooth[10] - f_smooth[0] ) / ( w_smooth[10] - w_smooth[0] )
#        w0_lin = array( [0.0, w_smooth[10] ], dtype = 'float_' )
#        f0_lin = array( [0.0, w_smooth[10] * secant_stiffness_w10 ], dtype = 'float_' )

        #axes.plot( w0_lin, f0_lin, color = 'black' )


    #--------------------------------------------------------------------------------
    # view
    #--------------------------------------------------------------------------------

    traits_view = View(VGroup(
                         Group(
                              Item('length', format_str = "%.3f"),
                              Item('width', format_str = "%.3f"),
                              Item('thickness', format_str = "%.3f"),
                              label = 'geometry'
                              ),
                         Group(
                              Item('loading_rate'),
                              Item('age'),
                              label = 'loading rate and age'
                              ),
                         Group(
                              Item('E_c', show_label = True, style = 'readonly', format_str = "%.0f"),
                              Item('ccs@', show_label = False),
                              label = 'composite cross section'
                              )
                         ),
                        scrollable = True,
                        resizable = True,
                        height = 0.8,
                        width = 0.6
                        )

if __name__ == '__main__':

    from matresdev.db.exdb.ex_run_table import ExRunClassExt
    ex = ExRunClassExt(klass = ExpBT3PT)
    ex.configure_traits()
