'''
Created on May 22, 2013

@author: rch
'''

import numpy as np
import pylab as p

offset = 0.1
arr = np.linspace(0 + offset, np.pi - offset)
print arr

tan_arr = 1.0 / np.tan(arr)

p.plot(arr, tan_arr)
p.show()
