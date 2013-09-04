'''
Created on Sep 2, 2013

@author: rch

The script uses the optimization methods to render the catenary
consisting of linear chain segments. The intention is to
implement the principle of minimum potential energy
numerically using the standard scipy methods.
'''

import numpy as np

# discretize the segments to get the initial position of the catenary

n_X = 4
x = np.linspace(0, 3, n_X)
y = np.zeros(n_X)
X = np.vstack([x, y]).T
L = np.array([[0, 1], [1, 2], [2, 3], [3, 4]])

print 'lengths', np.sqrt(np.sum((X[L[1:, 0]] - X[L[:-1, 0]]) ** 2, axis=1))

print 'xxxxxxx'
