'''
'''

import numpy as np
import sympy as sp

def index_mapping_node(dNx_mtx):    
    
    print 'dNx_mtx =\n', dNx_mtx
    dN_idx_map_one = ((0, 1, 1, 0), slice(None))
    print 'dNx_mtx, idx = \n', dNx_mtx[dN_idx_map_one]
     
    B_mtx = np.zeros ((3, 4, 2), dtype = 'f')
    B_enum_mtx = np.arange(24, dtype = 'f').reshape(3, 4, 2)
    B_idx_map_one = ((0, 1, 2, 2), slice(None), (0, 1, 0, 1))

    print 'B_mtx - orig\n', B_enum_mtx
    print 'B_mtx - selection\n', B_enum_mtx[B_idx_map_one]
    
    B_mtx[B_idx_map_one] = dNx_mtx[dN_idx_map_one]
    
    return B_mtx.reshape(3, 8)

''''Beispiel1'''

def index_mapping_four_Node (dNx_mtx):
    
    B_idx_map_four = ((0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2), (0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7))
     
    dN_idx_map_four = ((0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0), (0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3))
     
    B = np.zeros ((3, 8), dtype = 'f')
          
    B = np.zeros ((3, 8), dtype = 'f')

    B[B_idx_map_four] = dNx_mtx[dN_idx_map_four]
    
    return B


'''Beispiel 2'''

print "-----Beispiel four Node-------"

dNx_mtx4 = np.array([[11, 21, 31, 41], [12, 22, 32, 42]], dtype = 'f')

print "dNx_mtx = ", dNx_mtx4

B_mtx1 = index_mapping_node(dNx_mtx4)

print "B-Matrix 1 = \n", B_mtx1

print "B-Matrix 2 = \n", index_mapping_four_Node(dNx_mtx4)

print "-----Ende-----"

'''Ende Beispiel 2'''
