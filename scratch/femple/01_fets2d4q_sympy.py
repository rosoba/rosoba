
import numpy as np
import sympy as sp

x_, y_, z_, t_ = sp.symbols('x,y,z,t')

P = [1, x_, y_, x_ * y_]
PX = sp.lambdify([x_, y_], P)

# Nodal points
NP = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], dtype = 'f')

C = np.array([ PX(xi[0], xi[1]) for xi in NP ], dtype = 'f')

C1 = np.linalg.inv(C)

P_arr = np.array(P).reshape(1, 4)

N_mtx = sp.Matrix(np.dot(P_arr, C1))
N_mtx.simplify()

dN_mtx = N_mtx.jacobian([x_, y_]).T
dN_mtx.simplify()

N_fn = sp.lambdify([x_, y_], N_mtx)
dN_fn = sp.lambdify([x_, y_], dN_mtx)

print N_mtx
print dN_mtx
print N_fn(-1, -1)
print dN_fn(-1, -1)

# Now introduce the kinematics

# Use index operators to define the 

'''
'''

