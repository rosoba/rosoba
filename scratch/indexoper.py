'''
Created on Nov 25, 2013

@author: rch
'''

import numpy as np

L = np.array([[0, 1],
            [1, 2],
            [2, 3]], dtype='i')

X = np.array([[0, 0, 0],
             [1, 0, 0],
             [1, 1, 0],
             [0, 1, 0]], dtype='f')

u_i, u_j = X[L.T]

print u_j
