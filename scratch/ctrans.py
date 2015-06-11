'''
Created on Aug 17, 2014

@author: rch
'''

import numpy as np

vec_to_transform = np.array([0, 1, 0], dtype='f')

theta1 = np.pi / 4.
theta2 = np.pi - np.pi / 4.
theta3 = np.pi

def get_chi_mtx(theta):
    v1 = np.array([np.cos(theta), np.sin(theta), 0], dtype='f')
    z0 = np.array([0, 0, 1], dtype='f')
    v2 = np.cross(z0, v1)
    v3 = np.cross(v1, v2)
    base = np.array([v1 / np.linalg.norm(v1),
                     v2 / np.linalg.norm(v2),
                     v3 / np.linalg.norm(v3)], dtype='f')
    return base

chi1 = get_chi_mtx(theta1)
chi2 = get_chi_mtx(theta2)
chi3 = get_chi_mtx(theta3)
print 'chi1', chi1
print 'transformed', np.dot(chi1, vec_to_transform)

print np.dot(np.dot(chi1, chi2), chi3)

ct = np.array([], dtype='f')
