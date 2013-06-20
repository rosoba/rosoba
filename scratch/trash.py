'''
Created on May 22, 2013

@author: rch
'''

from traits.api import HasStrictTraits, Float, Property, cached_property
import math

class Circle(HasStrictTraits):

    radius = Float(1.0, auto_set=False, enter_set=True)

    area = Property(Float, depends_on='radius')
    @cached_property
    def _get_area(self):
        return self.radius ** 2 * math.pi

# all_traits_features.py --- Shows primary features of the Traits
#                            package

from traits.api import Delegate, HasTraits, Instance, \
                                 Int, Str
class Parent (HasTraits):

    # INITIALIZATION: last_name' is initialized to '':
    last_name = Str('')


class Child (HasTraits):

    age = Int

    # VALIDATION: 'father' must be a Parent instance:
    father = Instance(Parent)

    # DELEGATION: 'last_name' is delegated to father's 'last_name':
    last_name = Delegate('father')

    # NOTIFICATION: This method is called when 'age' changes:
    def _age_changed (self, old, new):
        print 'Age changed from %s to %s ' % (old, new)

# Set up the example:
joe = Parent()
joe.last_name = 'Johnson'
moe = Child()
moe.father = joe

# DELEGATION in action:
print "Moe's last name is %s " % moe.last_name
# Result:
# Moe's last name is Johnson

# NOTIFICATION in action
moe.age = 10
# Result:
# Age changed from 0 to 10

# VISUALIZATION: Displays a UI for editing moe's attributes
# (if a supported GUI toolkit is installed)
moe.configure_traits()


#if __name__ == '__main__':
#    c = Circle(radius=4.0)
#    print c.area
#    print c.area
#    c.radius = 5.0
#    print c.area
#
#    c.configure_traits()
