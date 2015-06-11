'''
Created on Feb 18, 2014

@author: tonda
'''

from traits.api import HasTraits, Button, Int
from traitsui.api import View, Item

class Calculator(HasTraits):

    enterx = Button

    number = Int

    view = View('enterx', 'number')

c1 = Calculator()

c1.configure_traits()
