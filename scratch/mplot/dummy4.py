from mayavi import mlab
import os.path

import pylab as pl


def mlab_imshowColor(im, alpha=255, **kwargs):
    """
    Plot a color image with mayavi.mlab.imshow.
    im is a ndarray with dim (n, m, 3) and scale (0->255]
    alpha is a single number or a ndarray with dim (n*m) and scale (0->255]
    **kwargs is passed onto mayavi.mlab.imshow(..., **kwargs)
    """
    try:
        alpha[0]
    except:
        alpha = pl.ones(im.shape[0] * im.shape[1]) * alpha
    if len(alpha.shape) != 1:
        alpha = alpha.flatten()

    # The lut is a Nx4 array, with the columns representing RGBA
    # (red, green, blue, alpha) coded with integers going from 0 to 255,
    # we create it by stacking all the pixles (r,g,b,alpha) as rows.
    myLut = pl.c_[im.reshape(-1, 3), alpha]
    myLutLookupArray = pl.arange(
        im.shape[0] * im.shape[1]).reshape(im.shape[0], im.shape[1])

    # We can display an color image by using mlab.imshow, a lut color list and
    # a lut lookup table.
    # temporary colormap
    theImshow = mlab.imshow(myLutLookupArray, colormap='binary', **kwargs)
    theImshow.module_manager.scalar_lut_manager.lut.table = myLut
    mlab.draw()

    return theImshow


def test_mlab_imshowColor():
    """
    Test if mlab_imshowColor displays correctly by plotting the wikipedia png example image
    """
    upath = os.path.expanduser('~')
    fname = os.path.join(
        upath, 'Downloads', 'eft01.png')
    im = pl.imread(fname, format='png')
    #im *= 255

    mlab_imshowColor(im[:, :, :3], im[:, :, -1])

if __name__ == "__main__":
    test_mlab_imshowColor()
