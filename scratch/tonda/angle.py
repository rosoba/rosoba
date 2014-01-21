'''
Created on Oct 28, 2013

@author: rch
'''

import math
import numpy as np

print 'pi =', math.pi, 'ty vole'

theta_1 = 2. / 3.*math.pi
theta_2 = 1. / 3.*math.pi
theta_3 = 1. / 3.*math.pi
theta_4 = 2. / 3.*math.pi

print 'suma theta =', theta_1 + theta_2 + theta_3 + theta_4 - 2 * math.pi, 'ty vole'

theta = [theta_1, theta_2, theta_3, theta_4, 1.4, -1.4 ]

print 'seznam theta', theta

for theta_i in theta:
    print 'theta', theta_i

print 'sum theta', np.sum(theta) - 2 * math.pi, 'ty vole'

