'''
Created on Jan 21, 2014

@author: rch
'''

from traits.api import HasStrictTraits, \
    Float, Property, cached_property, Int, \
    Trait, Event, on_trait_change, Instance, Button, Callable, \
    DelegatesTo, Constant, WeakRef, List

from util.traits.editors.mpl_figure_editor import \
    MPLFigureEditor

from matplotlib.figure import \
    Figure

from etsproxy.traits.ui.api import \
    View, Item, Group, HSplit, VGroup, HGroup

import numpy as np

from traitsui.api import \
    TableEditor, ObjectColumn, Label

bar_lst_editor = TableEditor(
                    columns=[ObjectColumn(name='x', label='x',
                                           editable=True,
                                           horizontal_alignment='center'),
                             ObjectColumn(name='z', label='z',
                                           editable=True,
                                           horizontal_alignment='center'),
                             ObjectColumn(name='A', label='A',
                                           editable=True,
                                           horizontal_alignment='center'),
                              ],
                    selection_mode='row',
                    selected='object.selected_var',
                    deletable=True,
                    editable=True,
                    show_toolbar=True,
                    auto_add=True,
                    configurable=True,
                    sortable=True,
                    reorderable=True,
                    sort_model=False,
                    orientation='vertical',
                    auto_size=True,
        )

class ECBCrossSectionGeo(HasStrictTraits):
    '''Base class for cross section types.
    '''

    cs_state = WeakRef

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

    n_layers = Int(10, auto_set=False, enter_set=True, geo_input=True)
    height = Float(0.3, auto_set=False, enter_set=True, geo_input=True)
    width = Float(0.2, auto_set=False, enter_set=True, geo_input=True)

    s_tex_z = Property(depends_on='+geo_input')
    '''property: spacing between the layers [m]
    '''
    @cached_property
    def _get_s_tex_z(self):
        return self.height / (self.n_layers + 1)

    z_ti_arr = Property(depends_on='+geo_input')
    '''property: distance from the top of each reinforcement layer [m]:
    '''
    @cached_property
    def _get_z_ti_arr(self):
        return np.array([ self.height - (i + 1) * self.s_tex_z
                         for i in range(self.n_layers) ],
                      dtype=float)

    zz_ti_arr = Property
    '''property: distance of reinforcement layers from the bottom
    '''
    def _get_zz_ti_arr(self):
        return self.height - self.z_ti_arr

    def _get_eps_ti_arr(self):
        return self.cs_state.eps_lo * self.zz_ti_arr

    def plot_geometry(self, ax):
        '''Plot geometry'''
        dx, dy = -self.width / 2, -self.height / 2
        xdata = np.array([0, self.width, self.width, 0, 0], dtype=float)
        ydata = np.array([0, 0, self.height, self.height, 0], dtype=float)

        ax.plot(xdata + dx, ydata + dy, color='blue')
        ax.hlines(self.zz_ti_arr + dy, [dx], [dx + self.width], lw=2, color='red')
        ax.axis([dx - 0.1 * self.width, dx + 1.1 * self.width, dy - 0.1 * self.height, dy + 1.1 * self.height])

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

class ECBRBar(HasStrictTraits):
    '''Class defining the position and cross sectional properties
    of a bar reinforcement.
    '''
    x = Float(0.0, auto_set=False, enter_set=True, geo_input=True)
    z = Float(0.0, auto_set=False, enter_set=True, geo_input=True)
    A = Float(0.0, auto_set=False, enter_set=True, geo_input=True)

    view = View(HGroup(
                Group(Item('x', springy=True),
                      Item('y'),
                      Item('A'),
                      label='Position and Area',
                      springy=True
                      ),
                springy=True,
                ),
                width=0.2,
                height=0.2,
                resizable=True,
                buttons=['OK', 'Cancel'])

class ECBGeoRB(ECBCrossSectionGeo):
    '''Rectangular cross section with arbitrarily spread bar reinforcement.
    '''
    #total height of crosssection
    #    
    height = Float(0.3, auto_set=False, enter_set=True, geo_input=True)

    #total width of crosssection
    #
    width = Float(0.2, auto_set=False, enter_set=True, geo_input=True)

    bar_coords = List(ECBRBar, [])
    '''List of reinforcement bars
    '''

    bar_coord_arr = Property
    def _get_bar_coord_arr(self):
        return np.array(self.bar_lst, dtype='f')

    #number of reinforcement bars
    #
    n_bars = Int(8, auto_set=False, enter_set=True, geo_input=True)

    #diameter of reinforcement bar
    #
    diameter_bar = Float(0.01, auto_set=False, enter_set=True, geo_input=True)

    A_bar = Property(depends_on='+geo_input')
    '''property: cross section area of one reinforcement bar [m**2]
    '''
    @cached_property
    def _get_A_bar(self):
        return (self.diameter_bar ** 2) / 4

    x_bar_arr = Property(depends_on='+geo_input')
    '''property: horizontal distance of individual reinforcement bars from left rim of cross section [m]
    '''
    def _get_x_bar_arr(self):
        return self.bar_coord_arr[:, 0]

    z_bar_arr = Property(depends_on='+geo_input')
    '''property: vertical distance of individual reinforcement bars from lower rim of cross section [m]
    '''
    def _get_z_bar_arr(self):
        return np.random.rand(self.n_bars) * self.height

    def plot_geometry(self, ax):
        '''Plot geometry'''
        '''Plot geometry'''
        dx, dy = -self.width / 2, -self.height / 2
        xdata = np.array([0, self.width, self.width, 0, 0], dtype=float)
        ydata = np.array([0, 0, self.height, self.height, 0], dtype=float)

        ax.plot(xdata + dx, ydata + dy, color='blue')
        ax.axis([dx - 0.1 * self.width, dx + 1.1 * self.width, dy - 0.1 * self.height, dy + 1.1 * self.height])

        ax.plot(self.x_bar_arr + dx, self.y_bar_arr + dy, 'ro')

    view = View(HSplit(Group(
                HGroup(
                Group(Item('width', springy=True),
                      Item('height'),
                      Item('@bar_lst', editor=bar_lst_editor),
                      Item('n_bars'),
                      Item('diameter_bar'),
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

class ECBGeoRI(ECBCrossSectionGeo):
    '''Rectangular cross section with arbitrarily spread bar reinforcement.
    '''

    #total height of crosssection
    #
    height = Float(0.6, auto_set=False, enter_set=True, geo_input=True)

    #height of upper flange
    #
    height_up = Float(0.1, auto_set=False, enter_set=True, geo_input=True)

    #height of lower flange
    #
    height_lo = Float(0.1, auto_set=False, enter_set=True, geo_input=True)

    #total width of crosssection
    #
    width = Float(0.6, auto_set=False, enter_set=True, geo_input=True)

    #width of stalk
    #
    width_st = Float(0.2, auto_set=False, enter_set=True, geo_input=True)

    #number of reinforcement bars
    #
    n_bars = Int(8, auto_set=False, enter_set=True, geo_input=True)

    #vertical distance of reinforcement bars from lower rim
    #
    z_bars = Float(0.05, auto_set=False, enter_set=True, geo_input=True)

    #diameter of reinforcement bar
    #
    diameter_bar = Float(0.01, auto_set=False, enter_set=True, geo_input=True)

    A_bar = Property(depends_on='+geo_input')
    '''property: cross section area of one reinforcement bar [m**2]
    '''
    @cached_property
    def _get_A_bar(self):
        return (self.diameter_bar ** 2) / 4

    x_bar_arr = Property(depends_on='+geo_input')
    '''property: horizontal distance of individual reinforcement bars from left rim of cross section [m]
    '''
    def _get_x_bar_arr(self):
        return np.array([ (i + 1) * self.width / (self.n_bars + 1)
                         for i in range(self.n_bars) ],
                      dtype=float)

    def plot_geometry(self, ax):
        '''Plot geometry'''
        dx, dy = -self.width / 2, -self.height / 2

        xdata = np.array([0, self.width, self.width, (self.width + self.width_st) / 2,
                          (self.width + self.width_st) / 2, self.width, self.width, 0, 0,
                          (self.width - self.width_st) / 2, (self.width - self.width_st) / 2,
                          0, 0], dtype=float)
        ydata = np.array([0, 0, self.height_lo, self.height_lo, self.height - self.height_up,
                          self.height - self.height_up, self.height, self.height, self.height - self.height_up,
                          self.height - self.height_up, self.height_lo, self.height_lo, 0], dtype=float)

        ax.plot(xdata + dx, ydata + dy, color='blue')

        #array of z coordinates of bars
        z_bar_arr = np.empty(self.n_bars); z_bar_arr.fill(self.z_bars)

        ax.plot(self.x_bar_arr + dx, z_bar_arr + dy, 'ro')
        ax.axis([dx - 0.1 * self.width, dx + 1.1 * self.width, dy - 0.1 * self.height, dy + 1.1 * self.height])

    view = View(HSplit(Group(
                HGroup(
                Group(Item('width', springy=True),
                      Item('width_st'),
                      Item('height'),
                      Item('height_lo'),
                      Item('height_up'),
                      Item('n_bars'),
                      Item('z_bars'),
                      Item('diameter_bar'),
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

class ECBGeoCB(ECBCrossSectionGeo):
    '''Circular cross section with bar reinforcement.
    '''

    #radius of cross section
    #
    radius = Float(0.3, auto_set=False, enter_set=True, geo_input=True)

    #number of reinforcement bars
    #
    n_bars = Int(10, auto_set=False, enter_set=True, geo_input=True)

    #distance of reinforcement bars from outer rim
    #
    d_bars = Float(0.2, auto_set=False, enter_set=True, geo_input=True)

    #diameter of reinforcement bar
    #
    diameter_bar = Float(0.01, auto_set=False, enter_set=True, geo_input=True)

    A_bar = Property(depends_on='+geo_input')
    '''property: cross section area of one reinforcement bar [m**2]
    '''
    @cached_property
    def _get_A_bar(self):
        return (self.diameter_bar ** 2) / 4

    def plot_geometry(self, ax):
        '''Plot geometry'''
        fi_outline_arr = np.append(np.arange(0, 2 * np.pi, np.pi / 60, dtype=float), 0.0)
        fi_bars_arr = np.arange(0, 2 * np.pi, 2 * np.pi / self.n_bars, dtype=float)

        ax.plot(np.cos(fi_outline_arr) * self.radius, np.sin(fi_outline_arr) * self.radius, color='blue')
        ax.plot(np.cos(fi_bars_arr) * (self.radius - self.d_bars),
                np.sin(fi_bars_arr) * (self.radius - self.d_bars), 'ro')
        ax.axis([-self.radius * 1.1, self.radius * 1.1, -self.radius * 1.1, self.radius * 1.1])

    view = View(HSplit(Group(
                HGroup(
                Group(Item('radius', springy=True),
                      Item('n_bars'),
                      Item('d_bars'),
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
    #ecs = ECBGeoRL(n_layers=4, height=0.2, width=0.1)

    #ecs.configure_traits()
    ecs = ECBGeoRB(height=0.3, width=0.2,
                   bar_coords=[[0.0, 0.6],
                               [0.2, 0.5]],
                   bar_areas=[0.04, 0.05, 0.04]
                   )

    ecs.configure_traits()
    #ecs = ECBGeoRI(height=0.6, width=0.6, height_lo = 0.1, height_up = 0.15, width_st = 0.08,
    #               n_bars = 8, z_bars = 0.05, diameter_bar = 0.01)

    #ecs.configure_traits()
    #ecs = ECBGeoCB(radius = 0.3, n_bars = 12, d_bars = 0.05)

    #ecs.configure_traits()
