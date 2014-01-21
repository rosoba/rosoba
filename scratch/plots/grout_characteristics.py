'''
Created on Nov 7, 2013

@author: rch
'''

import numpy as np
import mayavi.mlab as mlab

sig = np.array([[1.75, 1.88, 2.79],
                [1.56, 1.81, 2.13],
                [1.18, 4.11, 3.77]], dtype='f')
#sig_max = np.max(sig)
#sig /= sig_max

kappa = np.array([[1.4519, 2.5000, 6.7886],
                  [5.2363, 13.6182, 16.9232],
                  [14.4421, 22.9346, 39.0977]], dtype='f')
#kmax = np.max(kappa)
#kappa /= kmax

thickness = np.array([8, 14, 20], dtype='f')[:, None] * np.ones_like(sig)
angle = np.array([0, 10, 20], dtype='f')[None, :] * np.ones_like(sig)


mlab.figure(bgcolor=(1, 1, 1))
#mlab.barchart(thickness, angle, kappa, lateral_scale=0.4)
mlab.barchart(thickness, angle, sig, lateral_scale=0.4)

mlab.axes(xlabel='thickness',
          ylabel='angle',
          zlabel='crack initiation stress', nb_labels=3)

mlab.show()
