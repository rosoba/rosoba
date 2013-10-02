'''
Created on Sep 11, 2013

@author: rch

This sheet shows the evaluation of the potential energy
mgh of a discretization.
'''

import numpy as np

def print_where(arr, name='arr'):
    print name
    print 'positive'
    for i in np.where(arr == 1):
        print i
    print 'negative'
    for i in np.where(arr == -1):
        print i

#===============================================================================
# Index operator ... tensor calculus 
#===============================================================================

delta = np.zeros((3, 3,), dtype='f')
delta[(0, 1, 2), (0, 1, 2)] = 1

eps = np.zeros((3, 3, 3), dtype='f')
eps[(0, 1, 2), (1, 2, 0), (2, 0, 1)] = 1
eps[(2, 1, 0), (1, 0, 2), (0, 2, 1)] = -1

#===============================================================================
# Integration scheme
#===============================================================================

eta_ip = np.array([[1. / 3., 1. / 3.]], dtype='f')
eta_w = np.array([1. / 2.], dtype='f')

#===============================================================================
# Shape functions and their derivatives
#===============================================================================
def N(eta):
    return np.array([eta[:, 0], eta[:, 1], 1 - eta[:, 0] - eta[:, 1]], dtype='f').T

def N_deta(eta):
    return np.array([[[1, 0, -1],
                      [0, 1, -1]],
                     ], dtype='f')

N_eta_ip = N(eta_ip)
N_deta_ip = N_deta(eta_ip)

def get_E(x, F, debug=False):

    x_F = x[F]

    r = np.einsum('aK,IKi->Iai', N_eta_ip, x_F)
    r_deta = np.einsum('ajK,IKi->Iaij', N_deta_ip, x_F)
    n = np.einsum('Iai,Iaj,ijk->Iak', r_deta[..., 0], r_deta[..., 1], eps)
    a = np.sqrt(np.einsum('Iai,Iai->Ia', n, n))
    A = np.einsum('a,Ia->I', eta_w, a)
    E_I = np.einsum('a,Ia,Ia->I', eta_w, r[..., 2], a)
    E = np.sum(E_I)
    NN_delta_eps_x1 = np.einsum('aK,aL,KJ,jli,ILl->IaJji',
                                N_deta_ip[:, 0, :], N_deta_ip[:, 1, :], delta, eps, x_F)
    NN_delta_eps_x2 = np.einsum('aK,aL,LJ,kji,IKk->IaJji',
                                N_deta_ip[:, 0, :], N_deta_ip[:, 1, :], delta, eps, x_F)
    n_dx = NN_delta_eps_x1 + NN_delta_eps_x2
    a_dx = np.einsum('Ia,Iak,IaJjk->IaJj', 1 / a, n, n_dx)
    r3_a_dx = np.einsum('Ia,IaJj->IaJj', r[..., 2], a_dx)
    r3_dx = np.einsum('aK,KJ,j->aJj', N_eta_ip, delta, delta[2, :])
    a_r3_dx = np.einsum('Ia,aJj->IaJj', a, r3_dx)
    E_dx = np.einsum('a,IaJj->IJj', eta_w, (a_r3_dx + r3_a_dx))

    dof_map = (3 * F[:, :, np.newaxis] + np.arange(3)[np.newaxis, np.newaxis, :])

    E_dX = np.bincount(dof_map.flatten(), weights=E_dx.flatten())

    if debug:
        print
        print_where(eps, 'eps: Levi-Civita symbol indexes')
        print
        print 'N_eta_ip: shape function value in integration points\n', N_eta_ip
        print
        print 'N_deta_ip: local derivatives of shape function values in integration points\n', N_deta_ip
        print
        print 'r: global position of integration points\n', r
        print
        print 'r_deta: derivative of global ip position with respect to nodal coords\n', r_deta
        print
        print 'n: normal vectors in element I and integration point J', n
        print
        print 'a: area in element I and integration point J\n', a
        print
        print 'A: integrated area of element I\n', A
        print
        print 'E_I: potential energy of element I\n', E_I

        #print_where(NN_delta_eps_x1, 'NN_delta_eps_x1')
        #print_where(NN_delta_eps_x2, 'NN_delta_eps_x2')
        #print_where(n_dx, 'n_x')

        print
        print 'E: total potential energy of the crease pattern', E
        print
        print 'a_x: derivative of area with respect to nodal coordinates in element I integration point J \n', a_dx
        print
        print 'r3_a_dx\n', r3_a_dx
        print
        print 'r3_dx\n', r3_dx
        print
        print 'a_r3_dx\n', a_r3_dx
        print
        print 'E_dx\n', E_dx

    return E, E_dX

if __name__ == '__main__':
    #===============================================================================
    # Input nodes (x) and facets (F) 
    #===============================================================================
    x = np.array([[0, 0, 0],
                  [1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1],
                  [1, 0, 1],
                  [1, 1, 2],
                  [-1, 0, 0]], dtype='f')

    F = np.array([[0, 1, 2], [3, 4, 5], [1, 2, 6]], dtype='int')

    E, E_dx = get_E(x, F)
    print 'E', E
    print 'E_dx\n', E_dx


