from mayavi import mlab
import os.path

from oricreate.util.eulerangles import mat2euler, euler2mat
from oricreate.util.transformations import euler_matrix, euler_from_matrix
import pylab as pl

upath = os.path.expanduser('~')
fname = os.path.join(upath, 'Downloads', 'eft_01.png')
im = pl.imread(fname, format='png')[:, :, 0] * 255  # 1 color channel
rot = pl.r_[30, 80, 230]

R_orig = euler_matrix(*(rot * pl.pi / 180))[:3, :3]

#R_orig = euler2mat(*(rot * pl.pi / 180))
# print 'R_orig', R_orig
R_orig = pl.array([[0.89442719,  0.,    -0.4472136],
                   [-0.2981424,   0.74535599, -0.59628479],
                   [0.33333333,  0.66666667,  0.66666667]], dtype='float').T

print 'R_orig', R_orig
mlab.figure(bgcolor=(1, 1, 1), fgcolor=(0, 0, 0))
RList = [R_orig]  # , R_orig.T]
#
inorders = ['sxyz', 'sxzx', 'syxz', 'szxz', 'rzyx', 'rxzx', 'rzxy', 'rzxz', 'sxyx', 'syzx', 'syxy',
            'szyx', 'rxyx', 'rxzy', 'ryxy', 'rxyz', 'sxzy', 'syzy', 'szxy', 'szyz', 'ryzx', 'ryzy', 'ryxz', 'rzyz']

outorders = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]
for ii, inOrder in enumerate(inorders[2:3]):
    tries = 0
    # ]:
    for outOrder in outorders[2:3]:
        for R in RList:
            print 'R', R
            for vector, color in zip([[100, 0, 0], [0, 100, 0], [0, 0, 100]],
                                     [(1., 0., 0.), (0., 1., 0.), (0., 0., 1.)]):
                c = pl.c_[[0, tries * 1000, ii * 1000]]
                if ii == 0 and tries == 0:
                    vector = pl.r_[vector] * 5  # point out the first vector
                lin = R_orig.dot(pl.c_[[0, 0, 0], vector]) + c
                mlab.plot3d(*lin,
                            color=color,
                            tube_radius=5)

            lin3D = mlab.imshow(im, colormap="gray")
            print 'inOrder', inOrder
            rxyz = pl.array(
                euler_from_matrix(R, inOrder)) * 180 / pl.pi
            print 'outOrder', outOrder
            i, j, k = outOrder
            print 'o', [rxyz[i], rxyz[j], rxyz[k]]

            lin3D.actor.orientation = [rxyz[i], rxyz[j], rxyz[k]]
            lin3D.actor.position = c.flatten()
            lin3D.actor.scale = (6, 6, 1)
            tries += 1


mlab.draw()
mlab.show()
