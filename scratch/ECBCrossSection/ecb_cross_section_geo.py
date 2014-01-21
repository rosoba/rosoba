'''
Created on Jan 21, 2014

@author: rch
'''

from traits.api import HasStrictTraits

from etsproxy.traits.api import \
    HasStrictTraits, Float, Property, cached_property, Int, \
    Trait, Event, on_trait_change, Instance, Button, Callable, \
    DelegatesTo, Constant

from util.traits.editors.mpl_figure_editor import \
    MPLFigureEditor

from matplotlib.figure import \
    Figure

from etsproxy.traits.ui.api import \
    View, Item, Group, HSplit, VGroup, HGroup

import numpy as np

class ECBCrossSectionGeo(HasStrictTraits):
    '''Base class for cross section types.
    '''

    #===========================================================================
    # Plotting of the cross section
    #===========================================================================
    figure = Instance(Figure)
    def _figure_default(self):
        figure = Figure(facecolor='white')
        figure.add_axes([0.08, 0.13, 0.85, 0.74])
        return figure

    data_changed = Event

    replot = Button
    def _replot_fired(self):

        self.figure.clear()
        fig = self.figure
        ax = fig.add_subplot(111)

        self.plot_geometry(ax)

        self.data_changed = True

    def plot_geometry(self, ax):
        '''Plot geometry'''



class ECBGeoRL(ECBCrossSectionGeo):
    '''Rectangular cross section with layered fabric reinforcement.
    '''

    n_layers = Int(10, auto_set=False, enter_set=True)
    height = Float(0.3, auto_set=False, enter_set=True)
    width = Float(0.2, auto_set=False, enter_set=True)

    def plot_geometry(self, ax):
        '''Plot geometry'''
        xdata = [0, self.width, self.width, 0, 0]
        ydata = [0, 0, self.height, self.height, 0]
        ax.plot(xdata, ydata, color='blue')

    view = View(HSplit(Group(
                HGroup(
                Group(Item('n_layers', springy=True),
                      Item('height'),
                      Item('width'),
                      label='Geometry',
                      springy=True
                      ),
                springy=True,
                ),
                scrollable=True,
                ),
                Group(Item('replot', show_label=False),
                      Item('figure', editor=MPLFigureEditor(),
                           resizable=True, show_label=False),
                      id='simexdb.plot_sheet',
                      label='plot sheet',
                      dock='tab',
                      ),
                       ),
                width=0.8,
                height=0.7,
                resizable=True,
                buttons=['OK', 'Cancel'])

if __name__ == '__main__':
    ecs = ECBGeoRL(n_layers=4, height=0.2, width=0.1)

    ecs.configure_traits()
