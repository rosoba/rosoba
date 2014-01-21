'''
Created on Jan 15, 2014

@author: rch
'''

from traits.api import \
    HasStrictTraits, Float, Property, cached_property

from math import sqrt

class Point(HasStrictTraits):

    x = Float(0.0, auto_set=False, enter_set=True)
    y = Float(0.0, auto_set=False, enter_set=True)

    def move(self, dx, dy):
        self.x += dx
        self.y += dy

    distance_from_origin = Property(Float, depends_on='x,y')
    @cached_property
    def _get_distance_from_origin(self):
        print 'calculating distance'
        return sqrt(self.x ** 2 + self.y ** 2)

    def __str__(self):
        return 'x: %g, y: %g\n' % (self.x, self.y)

if __name__ == '__main__':
    p = Point()
    print 'point', p
    p.move(0.6, 4)
    print 'moved point, p', p

    print 'dist', p.distance_from_origin
    p.move(0.6, 4)
    print 'dist', p.distance_from_origin

