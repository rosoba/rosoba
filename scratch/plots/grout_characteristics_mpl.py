'''
Created on Nov 7, 2013

@author: rch
'''

from matplotlib import rc
#rc('font', **{'family':'sans-serif', 'sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font', **{'family':'serif', 'serif':['Palatino']})
rc('text', usetex=True)

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

row_names = ['0', '10', '20']
column_names = ['8', '14', '20']

sig = np.array([[1.75, 1.88, 2.79],
                [1.56, 1.81, 2.13],
                [1.18, 4.11, 3.77]], dtype='f')

kappa = np.array([[1.4519, 2.5000, 6.7886],
                  [5.2363, 13.6182, 16.9232],
                  [14.4421, 22.9346, 39.0977]], dtype='f')

thickness = np.array([8, 14, 20], dtype='f')[:, None] * np.ones_like(sig)
angle = np.array([0, 10, 20], dtype='f')[None, :] * np.ones_like(sig)

data = sig

lx = len(data[0])            # Work out matrix dimensions
ly = len(data[:, 0])
xpos = np.arange(0, lx, 1)    # Set up a mesh of positions
ypos = np.arange(0, ly, 1)
xpos, ypos = np.meshgrid(xpos + 0.25, ypos + 0.25)

xpos = xpos.flatten()   # Convert positions to 1D array
ypos = ypos.flatten()
zpos = np.zeros(lx * ly)

dx = 0.5 * np.ones_like(zpos)
dy = dx.copy()
dz = data.flatten()

print xpos
print ypos
print dz

ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b')

ticksx = np.arange(0.5, 5, 1)
plt.xticks(ticksx, column_names)

ticksy = np.arange(0.6, 7, 1)
plt.yticks(ticksy, row_names)

ax.set_xlabel('thickness [mm]')
ax.set_ylabel('angle [$^o$]')
ax.set_zlabel('$\sigma$ [N/mm$^2$]')

plt.show()
