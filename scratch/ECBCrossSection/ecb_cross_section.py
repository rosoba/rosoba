'''
Created on Jan 21, 2014

# @todo: consititutive law for compression zone and tensile zone.

'''

from traits.api import HasStrictTraits, \
    Float, Property, cached_property, Int, \
    Trait, Event, on_trait_change, Instance, Button, Callable, \
    DelegatesTo, Constant, Enum

from ecb_cross_section_geo import \
    ECBCrossSectionGeo, ECBGeoRL, ECBGeoRB, ECBGeoRI, ECBGeoCB

import numpy as np

from ecb_law import \
    ECBLBase, ECBLSteel

from constitutive_law import \
    ConstitutiveLawModelView

from cc_law import \
    CCLawBase, CCLawBlock, CCLawLinear, CCLawQuadratic, CCLawQuad


class ECBCrossSection(HasStrictTraits):
    '''Base class for cross section types.
    '''
    geo_type = Trait('RB', dict(RL=ECBGeoRL,
                               RB=ECBGeoRB,
                               RI=ECBGeoRI,
                               CB=ECBGeoCB))
    '''Type of cross section geometry
    RL - Rectangular with layers of reinforcement
    RB - Rectangular with bar reinforcement
    RI - I-shaped with bar reinforcement
    CB - Circular with bar reinforcement
    '''

    geo = Property(Instance(ECBCrossSectionGeo), depends_on='geo_type')
    '''Geometry of the cross section
    '''
    @cached_property
    def _get_geo(self):
        return self.geo_type_()

    cc_law_type = Trait('linear', dict(constant=CCLawBlock,
                                         linear=CCLawLinear,
                                         quadratic=CCLawQuadratic,
                                         quad=CCLawQuad),
                        cc_input=True)
    '''Selector of the concrete compression law type
    ['constant', 'linear', 'quadratic', 'quad']
    '''

    cc_law = Property(Instance(CCLawBase), depends_on='+cc_input')
    '''Compressive concrete law corresponding to cc_law_type'''
    @cached_property
    def _get_cc_law(self):
        return self.cc_law_type_(f_ck=50., eps_c_u=0.0033, cs=self)

    ecb_law = Property(Instance(ECBLBase), depends_on='+tt_input')
    '''Effective crack bridge law corresponding to ecb_law_type'''
    @cached_property
    def _get_ecb_law(self):
        return ECBLSteel()

    modified = Event
    @on_trait_change('geo.modified,+eps_input,+tt_input,+cc_input, cc_modified, tt_modified')
    def set_modified(self):
        self.modified = True

    eps_lo = Float(-0.0001, auto_set=False, enter_set=True, eps_input=True)
    eps_up = Float(-0.0001, auto_set=False, enter_set=True, eps_input=True)

    unit_conversion_factor = Constant(1000.0)
    '''Convert the MN to kN
    '''

    x = Property(depends_on='+eps_input,geo.modified')
    '''Height of the compressive zone
    '''

    def _get_x(self):
        if self.eps_up == self.eps_lo:
            # @todo: explain
            return (self.geo.height)
        else:
            return (abs(self.eps_up) / (abs(self.eps_up - self.eps_lo)) *
                     self.geo.height)

    n_cj = Int(30, auto_set=False, enter_set=True,
                 cc_input=True, eps_input=True)
    '''Number of subdivisions of the compressive zone
    '''

    z_cj_arr = Property(depends_on='+eps_input,geo.modified')
    '''Discretizaton of the  compressive zone
    '''
    @cached_property
    def _get_z_cj_arr(self):
        if self.eps_up <= 0: # bending
            zx = min(self.geo.height, self.x)
            return np.linspace(0, zx, self.n_cj)
        elif self.eps_lo <= 0: # bending
            return np.linspace(self.x, self.geo.height, self.n_cj)
        else: # no compression
            return np.array([0], dtype='f')

    w_cj_arr = Property(depends_on='geo.modified,z_cj_arr')
    '''Discretizaton of the  compressive zone - weight factors for general cross section
    '''
    @cached_property
    def _get_w_cj_arr(self):
        return self.geo.get_width(self.z_cj_arr)

    eps_cj_arr = Property(depends_on='+eps_input,geo.modified')
    '''Compressive strain at each integration layer of the compressive zone [-]:
    '''
    @cached_property
    def _get_eps_cj_arr(self):
        eps_j_arr = (self.eps_up + (self.eps_lo - self.eps_up) * self.z_cj_arr /
                     self.geo.height)
        return (-np.fabs(eps_j_arr) + eps_j_arr) / 2.0

    sig_cj_arr = Property(depends_on='+eps_input, +cc_input, cc_modified')
    '''Stresses at the j-th integration point.
    '''
    @cached_property
    def _get_sig_cj_arr(self):
        return -self.cc_law.mfn_vct(-self.eps_cj_arr)

    f_cj_arr = Property(depends_on='+eps_input, +cc_input, cc_modified, geo.modified')
    '''Layer force corresponding to the j-th integration point [kN].
    '''
    @cached_property
    def _get_f_cj_arr(self):
        return self.w_cj_arr * self.sig_cj_arr * self.unit_conversion_factor

    eps_bar_arr = Property(depends_on='+eps_input,geo.modified')
    '''Strain at the level of the i-th reinforcement bar
    '''
    @cached_property
    def _get_eps_bar_arr(self):
        return self.eps_lo + (self.eps_up - self.eps_lo) * self.geo.z_bar_arr / self.geo.height

    sig_bar_arr = Property(depends_on='+eps_input, geo.modified, +tt_input, tt_modified')
    '''Stresses at the i-th fabric layer.
    '''
    @cached_property
    def _get_sig_bar_arr(self):
        return self.ecb_law.mfn_vct(np.fabs(self.eps_bar_arr)) * np.sign(self.eps_bar_arr)

    f_bar_arr = Property(depends_on='+eps_input, geo.modified, +tt_input, tt_modified')
    '''force at the height of each reinforcement layer [kN]:
    '''
    @cached_property
    def _get_f_bar_arr(self):
        sig_bar_arr = self.sig_bar_arr
        A_bar_arr = self.geo.bar_area_arr
        return sig_bar_arr * A_bar_arr * self.unit_conversion_factor

    N = Property(depends_on='geo.modified,+eps_input,+tt_input,+cc_input, cc_modified, tt_modified')
    '''Get the resulting normal force.
    '''
    @cached_property
    def _get_N(self):
        N_tk = sum(self.f_bar_arr)
        N_ck = np.trapz(self.f_cj_arr, self.z_cj_arr)

        N_internal = N_ck + N_tk
        return N_internal

    M = Property(depends_on='geo.modified,+eps_input,+tt_input, +cc_input, cc_modified, tt_modified')
    '''Get the resulting moment.
    '''
    @cached_property
    def _get_M(self):
        M_tk = np.dot(self.f_bar_arr, (self.geo.height - self.geo.z_bar_arr))
        M_ck = np.trapz(self.f_cj_arr * self.z_cj_arr, self.z_cj_arr)
        M_internal_ = M_tk + M_ck

        M_internal = M_internal_ - self.N * self.geo.gravity_centre
        '''Moment evaluated with respect to the gravity centre
        '''
        return M_internal

if __name__ == '__main__':
    ecs = ECBCrossSection()

    ecs.geo_type = 'RB'
    print 'Moment for rectangular cross section with bar reinforcement', ecs.M
    print 'Normal force for rectangular cross section with bar reinforcement', ecs.N
    ecs.geo_type = 'RI'
    print 'Moment for I-shaped cross section with bar reinforcement', ecs.M
    print 'Normal force for I-shaped cross section with bar reinforcement', ecs.N
    ecs.geo_type = 'CB'
    print 'Moment for circular cross section with bar reinforcement', ecs.M
    print 'Normal force for circular cross section with bar reinforcement', ecs.N

    ecs.configure_traits()

