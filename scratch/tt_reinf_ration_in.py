'''
Created on Jun 9, 2015

@author: rch
'''
import numpy as np

E_f = 240000.0
V_f = np.linspace(0.001, 1.0, 11)
E_m = 28000.0
eps_c = 0.0001

E_c = E_f * V_f + E_m * (1 - V_f)
sig_f = (E_c - E_m * (1 - V_f)) / V_f * eps_c

print E_c
print sig_f
