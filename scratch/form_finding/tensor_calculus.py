'''
Created on Sep 11, 2013

@author: rch
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
print_where(eps, 'eps')


eta_ip = np.array([1. / 3., 1. / 3.], dtype='f')
eta_w = np.array([1. / 2.], dtype='f')

def N(eta):
    return np.array([eta[0], eta[1], 1 - eta[0] - eta[1]], dtype='f')

def N_deta(eta):
    return np.array([[1, 0, -1],
                     [0, 1, -1]], dtype='f')

N_eta_ip = N(eta_ip)
N_deta_ip = N_deta(eta_ip)
print 'N_eta_ip\n', N_eta_ip
print 'N_deta_ip\n', N_deta_ip

x = np.array([[[0, 0, 0], [1, 0, 0], [0, 1, 0]],
              [[0, 0, 1], [1, 0, 1], [1, 1, 2]],
              [[1, 0, 0], [0, 1, 0], [-1, 0, 0]],
              ], dtype='f')

r = np.einsum('k,gki->gi', N_eta_ip, x)
print 'r\n', r

r_deta = np.einsum('jk,gki->gij', N_deta_ip, x)
print 'r_deta\n', r_deta

n = np.einsum('gi,gj,ijk->gk', r_deta[..., 0], r_deta[..., 1], eps)
print 'n\n', n

a2 = np.einsum('gi,gi->g', n, n)
print 'a2\n', a2

a = np.sqrt(a2)
print 'a\n', a

A = eta_w * a
print 'A\n', A

E = eta_w * r[:, 2] * a
print 'E\n', E

delta_delta_eps = np.einsum('mi,oj,opk->mijpk', delta, delta, eps)
print_where(delta_delta_eps, 'delta_delta_eps')

delta_eps = np.einsum('mi,jpk->mijpk', delta, eps)
print_where(delta_eps, 'delta_eps')

delta_eps_x = np.einsum('mi,jpk,gnp->gmijkn', delta, eps, x)
print_where(delta_eps_x, 'delta_eps_x')

NN_delta_eps_x1 = np.einsum('m,n,mi,jpk,gnp->gijk', N_deta_ip[0, :], N_deta_ip[1, :], delta, eps, x)
print_where(NN_delta_eps_x1, 'NN_delta_eps_x1')

NN_delta_eps_x2 = np.einsum('m,n,ni,ojk,gmo->gijk', N_deta_ip[0, :], N_deta_ip[1, :], delta, eps, x)
print_where(NN_delta_eps_x2, 'NN_delta_eps_x2')

n_dx = NN_delta_eps_x1 + NN_delta_eps_x2
print_where(n_dx, 'n_x')

a_dx = np.einsum('g,gk,gijk->gij', 1 / a, n, n_dx)
print 'a_x\n', a_dx

r3_a_dx = np.einsum('g,gij->gij', r[:, 2], a_dx)
print 'r3_a_dx\n', r3_a_dx

r3_dx = np.einsum('k,ki,j->ij', N_eta_ip, delta, delta[2, :])
print 'r3_dx', r3_dx

a_r3_dx = np.einsum('g,gij->gij', a, r3_dx[np.newaxis, :])
print 'a_r3_dx', a_r3_dx

E_dx = eta_w * (a_r3_dx + r3_a_dx)
print 'E_dx', E_dx

