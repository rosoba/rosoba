'''
Created on Jun 17, 2014

@author: rch
'''

from scipy.stats.distributions import gamma
import numpy as np

x = np.linspace(0, 0.5, 4000)

shapes = np.array([0.18, 0.15, 0.23, 0.15, 0.20, 0.18], dtype='f')
scales = np.array([0.8, 0.7, 0.6, 1.0, 0.7, 0.8], dtype='f')
locs = np.array([0.0055, 0.0080, 0.0010, 0.0050, 0.0090, 0.0057], dtype='f')

import pylab as p
for i in range(0, len(shapes)):
    g = gamma(shapes[i], scale=scales[i], loc=locs[i])
    pdf = g.pdf(x)
    p.plot(x, pdf)

p.show()

