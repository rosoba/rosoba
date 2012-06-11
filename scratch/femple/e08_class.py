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

import numpy as np

from scipy.linalg import \
     inv

from ibvpy.fets.fets_eval import FETSEval
import math

  
class geotransform (HasTraits):

    debug_on = True
       '''
    L1 = 0.2
    L2 = 0.2
    alpha = math.pi / 2.0 / 3.0
    d = 0.01
    h = 0.01

    n_z = 4
    n_x = 5 
 
    def __init__(self,X,alpha,L2,d):
        self.X = X
        self.alpha = alpha
        self.L2=L2
        self.d=d
   '''
    gt = Property
    
    def _get_gt(self,L2, d, alpha, points):
        
        X1 = np.array([[L2, 0], [0, 0], [0, d], [-L2, d]], dtype = 'f')
        T = np.array([[ math.cos(alpha), math.sin(alpha)],
                  [ -math.sin(alpha), math.cos(alpha)]], dtype = 'f')
        X1 = np.dot(X1, T)
        return self.N_transform(points, X1)

    
    def N_transform(self,r, X):
        
        geo_r = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], dtype = 'f')
        cx = np.array(geo_r, dtype = 'float_')
        Nr = np.array([1 / 4. * (1 + r[:, 0] * cx[i, 0]) * (1 + r[:, 1] * cx[i, 1])
                      for i in range(0, 4) ])
        return np.dot(Nr.T, X)      