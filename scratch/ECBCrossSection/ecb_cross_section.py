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


class ECBCrossSection(HasStrictTraits):
    '''Base class for cross section types.
    '''
    geo_type = Trait('RL', dict(RL=ECBGeoRL,
                               RB=ECBGeoRB,
                               RI=ECBGeoRI,
                               CB=ECBGeoCB))
    '''Type of cross section geometry
    RL - Rectangular with layers of reinforcement
    '''

    geo = Property(Instance(ECBCrossSectionGeo), depends_on='geo_type')
    '''Geometry of the cross section
    '''
    @cached_property
    def _get_geo(self):
        return self.geo_type_(cs_state=self)

    eps_lo = Float()

if __name__ == '__main__':
    ecs = ECBCrossSection()

    print 'Moment for rectangular cross-section with layered reinforcement', ecs.M
    ecs.geo_type = 'CB'
    print 'Moment for circular cross section with bar reinforcement', ecs.M

    ecs.configure_traits()

