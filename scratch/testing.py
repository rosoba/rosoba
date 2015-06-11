'''
Created on Sep 25, 2014

@author: rch
'''


import numpy as np

def make_something_with_array():

    a = np.linspace(0, 20, 20)
    a_idx = np.arange(0, 20, 2, dtype='int')

    b_idx = a_idx.reshape(2, 5)

    return a[b_idx]


print make_something_with_array()
