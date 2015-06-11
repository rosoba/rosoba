'''
Created on Oct 28, 2013

@author: rch
'''

import math
import numpy as np
import pylab as p
#
#x = np.linspace(0.0001, 10, 100)
#
#print x
#
#y = np.log(x, 10)
#
#print y
#
#p.plot(y, x, color='red')
#p.plot(x, y, color='blue')
#
#p.show()

jacke_vorher = 229.0
jacke_jetzt = 200.0

procent_vom_urspruenlichen_preis = jacke_jetzt / jacke_vorher

procent_reduktion = 1 - procent_vom_urspruenlichen_preis

print procent_reduktion

import sympy

x_, jv_, jj_ = sympy.symbols('x, jv, jj')

gleichung = 100.0 / jv_ - x_ / jj_

print gleichung

print sympy.solve(gleichung, 'x')



