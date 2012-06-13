'''
Created on 08.06.2012

@author: demian
'''

import sys
print sys.path

from etsproxy.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property, \
     Constant

import etsproxy.traits.has_traits
etsproxy.traits.has_traits.CHECK_INTERFACES = 2

from etsproxy.traits.ui.api import \
     Item, View, HGroup, ListEditor, VGroup, Group

from etsproxy.traits.ui.menu import \
     NoButtons, OKButton, CancelButton, Action, CloseAction, Menu, \
     MenuBar, Separator

from ibvpy.api import IBVModel

import numpy as np

from scipy.linalg import \
     inv

from ibvpy.fets.fets_eval import FETSEval
import math

class FoldedBondTest(IBVModel):
    '''Idealization of the test for the characterization of
    bond behavior within crease line of a folded plate. 
    '''
    L1 = Float(0.2, desc = 'Length of the left part')
    L2 = Float(0.2, desc = 'Length of the right part')
    alpha = Float(math.pi / 2.0 / 3.0, desc = 'Fold angle')
    d = Float(0.01, desc = 'thickness of the plate')
    h = Float(0.01, desc = 'width of the plate')

    #===========================================================================
    # Discretization parameters
    #===========================================================================
    n_z = Int(4, desc = 'number of elements in the thickness direction')
    n_x = Int(5, desc = 'number of elements in the length direction of a plate')

    view = View(Item('L1', label = 'length of fixed part'),
                Item('L2', label = 'length of the loaded part'),
                Item('alpha'),
                width = 0.2,
                height = 0.3
                )

if __name__ == '__main__':
    fbt = FoldedBondTest()
    fbt.configure_traits()
