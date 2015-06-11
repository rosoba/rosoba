'''
Created on Apr 10, 2014

@author: rch


finite element code consists of the following components

array of nodes with [n_el, n_d]

'''

from ibvpy.mesh.fe_grid import \
    FEGrid

from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q

fets_sample = FETS2D4Q()
fe_domain = FEGrid(shape=(2, 3),
                   coord_max=(2, 3),
                   fets_eval=fets_sample)

import sys
print 'refcount', sys.getrefcount(fe_domain)
dots = fe_domain.dots
print 'dots', dots
print 'fets', dots.fets_eval
print 'refcount', sys.getrefcount(fe_domain)

print 'dof_r'
print fe_domain.dof_r

print fe_domain.geo_grid.cell_node_map
print fe_domain.dof_grid.cell_dof_map
print fe_domain.elem_dof_map
print fe_domain.elem_X_map[0]

#    for e in fe_domain.elements:
#        print e


