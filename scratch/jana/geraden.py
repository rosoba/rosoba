'''
Created on Dec 2, 2013

@author: rch
'''

import numpy as np
import sympy as sp

a, b, x, y = sp.symbols('a,b,x,y')

F = a * x + b - y

F.subs({x:0, y:3})
print F
F.subs({x:2, y:5})

print F

