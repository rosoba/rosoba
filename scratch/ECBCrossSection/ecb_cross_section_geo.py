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
    '''Number of textile reinforcement layers
    '''
    
    height = Float(0.3, auto_set=False, enter_set=True, geo_input=True)
    '''Height of cross section
    '''

    width = Float(0.2, auto_set=False, enter_set=True, geo_input=True)
    '''Width of cross section
    '''
    
    s_tex_z = Property(depends_on='+geo_input')
    '''property: spacing between the layers [m]
    '''
    @cached_property
    def _get_s_tex_z(self):
        return self.height / (self.n_layers + 1)

    z_ti_arr = Property(depends_on='+geo_input')
    '''property: distance of each reinforcement layer from the top [m]:
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

    def plot_geometry(self, ax):
        '''Plot geometry'''
        dx, dz = -self.width / 2, -self.height / 2
        xdata = np.array([0, self.width, self.width, 0, 0], dtype=float)
        zdata = np.array([0, 0, self.height, self.height, 0], dtype=float)

        ax.plot(xdata + dx, zdata + dz, color='blue')
        ax.hlines(self.zz_ti_arr + dz, [dx], [dx + self.width], lw=2, color='red')
        ax.axis([dx - 0.1 * self.width, dx + 1.1 * self.width, dz - 0.1 * self.height, dz + 1.1 * self.height])

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

    height = Float(0.3, auto_set=False, enter_set=True, geo_input=True)
    '''total height of crosssection
    '''

    width = Float(0.2, auto_set=False, enter_set=True, geo_input=True)
    '''total width of crosssection
    '''

    bar_coords = List([[0.05, 0.05], [0.15, 0.05], [0.05, 0.25], [0.15, 0.25]])
    '''list containing coordinates of individual reinforcement bars
    '''

    bar_areas = List([0.0004, 0.0004, 0.0004, 0.0004])
    '''list containing areas of individual reinforcement bars
    '''
    
    bar_lst = List(ECBRBar, [])
    '''List of reinforcement bars
    '''

    gravity_centre = Property(depends_on='geo_input')
    '''z distance of gravity centre from upper rim
    '''
    def _get_gravity_centre(self):
        return self.height / 2

    def get_width(self, z):
        width_arr = np.empty(z.shape)
        width_arr[:] = self.width
        return width_arr
    
    bar_area_arr = Property(depends_on='+geo_input')
    '''property: area of individual reinforcement bars [m**2]
    '''
    def _get_bar_area_arr(self):
        return np.array(self.bar_areas, dtype='f')
    
    bar_coord_arr = Property(depends_on='+geo_input')
    '''property: coordinates of individual reinforcement bars with respect to the left lower corner of the cross section [m]
    '''
    def _get_bar_coord_arr(self):
        return np.array(self.bar_coords, dtype='f')

    x_bar_arr = Property(depends_on='+geo_input')
    '''property: horizontal distance of individual reinforcement bars from left rim of cross section [m]
    '''
    def _get_x_bar_arr(self):
        return self.bar_coord_arr[:, 0]

    z_bar_arr = Property(depends_on='+geo_input')
    '''property: vertical distance of individual reinforcement bars from lower rim of cross section [m]
    '''
    def _get_z_bar_arr(self):
        return self.bar_coord_arr[:, 1]

    def plot_geometry(self, ax):
        '''Plot geometry'''
        dx, dz = -self.width / 2, -self.height / 2
        xdata = np.array([0, self.width, self.width, 0, 0], dtype=float)
        zdata = np.array([0, 0, self.height, self.height, 0], dtype=float)
        ax.plot(xdata + dx, zdata + dz, color='blue')
        ax.axis([dx - 0.1 * self.width, dx + 1.1 * self.width, dz - 0.1 * self.height, dz + 1.1 * self.height])

        ax.plot(self.x_bar_arr + dx, self.z_bar_arr + dz, 'ro')

    view = View(HSplit(Group(
                HGroup(
                Group(Item('width', springy=True),
                      Item('height'),
#                      Item('@bar_lst', editor=bar_lst_editor),
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

    height = Float(0.6, auto_set=False, enter_set=True, geo_input=True)
    '''total height of crosssection
    '''

    height_up = Float(0.1, auto_set=False, enter_set=True, geo_input=True)
    '''height of upper flange
    '''

    height_lo = Float(0.1, auto_set=False, enter_set=True, geo_input=True)
    '''height of lower flange
    '''

    width_up = Float(0.4, auto_set=False, enter_set=True, geo_input=True)
    '''width of upper flange
    '''

    width_lo = Float(0.6, auto_set=False, enter_set=True, geo_input=True)
    '''width of lower flange
    '''

    width_st = Float(0.2, auto_set=False, enter_set=True, geo_input=True)
    '''width of stalk
    '''
    
    bar_coords = List([[0.2, 0.05], [0.4, 0.05], [0.2, 0.55], [0.4, 0.55]])
    '''list containing coordinates of individual reinforcement bars
    '''

    bar_areas = List([0.0004, 0.0004, 0.0004, 0.0004])
    '''list containing areas of individual reinforcement bars
    '''

    gravity_centre = Property(depends_on='geo_input')
    '''z distance of gravity centre from upper rim
    '''
    def _get_gravity_centre(self):
        A_up, z_up = self.width_up * self.height_up, self.height_up / 2
        A_lo, z_lo = self.width_lo * self.height_lo, self.height - self.height_lo / 2
        A_st, z_st = self.width_st * (self.height - self.height_up - self.height_lo), (self.height + self.height_up - self.height_lo) / 2
        return (A_up*z_up + A_lo*z_lo + A_st*z_st) / (A_up + A_lo + A_st)
    
    def get_width(self, z):
        '''returns width of cross section for different vertical coordinates
        '''
        width_arr = self.width_up + (np.sign(z-self.height_up)+1)/2 * (self.width_st - self.width_up) + \
        (np.sign(z-self.height + self.height_lo)+1)/2 * (self.width_lo - self.width_st)
        return width_arr
        
    bar_area_arr = Property(depends_on='+geo_input')
    '''property: area of individual reinforcement bars [m**2]
    '''
    def _get_bar_area_arr(self):
        return np.array(self.bar_areas, dtype='f')

    bar_coord_arr = Property(depends_on='+geo_input')
    '''property: coordinates of individual reinforcement bars with respect to the left lower corner of the cross section [m]
    '''
    def _get_bar_coord_arr(self):
        return np.array(self.bar_coords, dtype='f')

    x_bar_arr = Property(depends_on='+geo_input')
    '''property: horizontal distance of individual reinforcement bars from left rim of cross section [m]
    '''
    def _get_x_bar_arr(self):
        return self.bar_coord_arr[:, 0]

    z_bar_arr = Property(depends_on='+geo_input')
    '''property: vertical distance of individual reinforcement bars from lower rim of cross section [m]
    '''
    def _get_z_bar_arr(self):
        return self.bar_coord_arr[:, 1]

    def plot_geometry(self, ax):
        '''Plot geometry'''
        w_max = max(self.width_lo, self.width_up)
       
        dx, dz = -w_max / 2, -self.height / 2

        xdata = np.array([- self.width_lo / 2, self.width_lo / 2, self.width_lo / 2, self.width_st / 2,
                          self.width_st / 2, self.width_up / 2, self.width_up / 2, - self.width_up / 2,
                          - self.width_up / 2, - self.width_st / 2, - self.width_st / 2,
                          - self.width_lo / 2, - self.width_lo / 2], dtype=float) + w_max / 2
                          
        zdata = np.array([0, 0, self.height_lo, self.height_lo, self.height - self.height_up,
                          self.height - self.height_up, self.height, self.height, self.height - self.height_up,
                          self.height - self.height_up, self.height_lo, self.height_lo, 0], dtype=float)

        ax.plot(xdata + dx, zdata + dz, color='blue')

        ax.plot(self.x_bar_arr + dx, self.z_bar_arr + dz, 'ro')
        ax.axis([dx - 0.1 * w_max, dx + 1.1 * w_max, dz - 0.1 * self.height, dz + 1.1 * self.height])

    view = View(HSplit(Group(
                HGroup(
                Group(Item('width_up', springy=True),
                      Item('width_lo'),
                      Item('width_st'),
                      Item('height'),
                      Item('height_lo'),
                      Item('height_up'),
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

    radius = Float(0.3, auto_set=False, enter_set=True, geo_input=True)
    '''radius of cross section
    '''

    n_bars = Int(10, auto_set=False, enter_set=True, geo_input=True)
    '''number of reinforcement bars
    '''

    d_bars = Float(0.2, auto_set=False, enter_set=True, geo_input=True)
    '''distance of reinforcement bars from outer rim
    '''

    A_bar = Float(0.0004, auto_set=False, enter_set=True, geo_input=True)
    '''diameter of reinforcement bar
    '''
    
    height = Property(depends_on='radius')
    '''Height of cross section
    '''
    def _get_height(self):
        return 2 * self.radius

    gravity_centre = Property(depends_on='geo_input')
    '''z distance of gravity centre from upper rim
    '''
    def _get_gravity_centre(self):
        return self.radius

    def get_width(self, z):
        r_dist_arr = z - self.radius
        '''distance of the point from center
        '''
        width_arr = 2 * np.sqrt(self.radius ** 2 - r_dist_arr ** 2)
        return width_arr
        
    fi_bar_arr = Property(depends_on='+geo_input')
    '''property: array containing degrees for computation of bar coordinates
    '''
    def _get_fi_bar_arr(self):
        return np.arange(0, 2 * np.pi, 2 * np.pi / self.n_bars, dtype=float)
    
    bar_area_arr = Property(depends_on='+geo_input')
    '''property: area of individual reinforcement bars [m**2]
    '''
    def _get_bar_area_arr(self):
        return np.ones(self.n_bars) * self.A_bar

    x_bar_arr = Property(depends_on='+geo_input')
    '''property: horizontal distance of individual reinforcement bars from left rim of cross section [m]
    '''
    def _get_x_bar_arr(self):
        return np.cos(self.fi_bar_arr) * (self.radius - self.d_bars) + self.radius

    z_bar_arr = Property(depends_on='+geo_input')
    '''property: vertical distance of individual reinforcement bars from lower rim of cross section [m]
    '''
    def _get_z_bar_arr(self):
        return np.sin(self.fi_bar_arr) * (self.radius - self.d_bars) + self.radius

    def plot_geometry(self, ax):
        '''Plot geometry'''
        fi_outline_arr = np.append(np.arange(0, 2 * np.pi, np.pi / 60, dtype=float), 0.0)

        ax.plot(np.cos(fi_outline_arr) * self.radius, np.sin(fi_outline_arr) * self.radius, color='blue')
        ax.plot(self.x_bar_arr - self.radius, self.z_bar_arr - self.radius, 'ro')
        ax.axis([-self.radius * 1.1, self.radius * 1.1, -self.radius * 1.1, self.radius * 1.1])

    view = View(HSplit(Group(
                HGroup(
                Group(Item('radius', springy=True),
                      Item('n_bars'),
                      Item('d_bars'),
                      Item('A_bar'),
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

    ecs = ECBGeoRI(width_up=0.8, width_lo=0.5)

    ecs.configure_traits()
