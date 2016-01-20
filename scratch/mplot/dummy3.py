
import mayavi.mlab as m
import numpy as np


def test_plot3d():
    """Generates a pretty set of lines."""
    n_mer, n_long = 6, 11
    pi = np.pi
    dphi = pi / 1000.0
    phi = np.arange(0.0, 2 * pi + 0.5 * dphi, dphi)
    mu = phi * n_mer
    x = np.cos(mu) * (1 + np.cos(n_long * mu / n_mer) * 0.5)
    y = np.sin(mu) * (1 + np.cos(n_long * mu / n_mer) * 0.5)
    z = np.sin(n_long * mu / n_mer) * 0.5

    l = m.plot3d(x, y, z, np.sin(mu), tube_radius=0.025, colormap='Spectral')

    return l

if __name__ == '__main__':

    m.text3d(1, 0, 0, 'function', orient_to_camera=False,
             orientation=(0, 0, 90), scale=0.4)

    m.show_pipeline()

    m.pipeline

    #    test_plot3d()
    m.show()
