'''
Created on Jan 13, 2014

@author: rch
'''

from traits.api import \
    HasStrictTraits, Int, Float, Instance

from traitsui.api import \
    View, Item, ModelView

from stats.misc.random_field.random_field_1D import \
    RandomField

class TensileSpecimen(HasStrictTraits):

    strength_field = Instance(RandomField)
    def _strength_field_default(self):
        return RandomField()

if __name__ == '__main__':

    sf = RandomField(length=0.6, lacor=0.05, distr_type='Weibull',
                     scale=4.4, stdev=0.5)
    ts = TensileSpecimen(strength_field=sf)
    rf = ts.strength_field.random_field
    x = ts.strength_field.xgrid

    import pylab as plt
    import numpy as np

    plt.plot(x, rf)
    plt.ylim(0, np.max(rf) * 1.1)
    plt.show()
