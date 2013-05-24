'''
Created on May 24, 2013

@author: rch
'''

'''
Example of a 2d arch
In section 'calculation parameters' the following predefined setting are available:
- loading type: symmetric (single displacement at top) applied at top
                asymmetric (single displacement at the left side of the arch applied radially at L/4;
- geometry: the opening angle of the arch can be specified (=< Pi);
            a varying thickness from bottom to top can be defined
            with either linear or quadratic transition;
- fets:     choose between linear of quadratic 2d elements: FETS2D4Q or FETS2D9Q;
'''

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
    TLine, BCDof, IBVPSolve as IS, DOTSEval
from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
from ibvpy.api import BCDofGroup
from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
from ibvpy.fets.fets2D.fets2D9q import FETS2D9Q
from ibvpy.mesh.fe_grid import FEGrid

from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import MATS2DScalarDamage
from ibvpy.mats.mats2D.mats2D_sdamage.strain_norm2d import Energy, Euclidean, Mises, Rankine, Mazars
from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm import \
    MATS2DMicroplaneDamage, PhiFnGeneral, PhiFnStrainHardening

from ibvpy.mats.mats_proxy import MATSProxy
from mathkit.geo.geo_ndgrid import GeoNDGrid
from mathkit.mfn.mfn_ndgrid.mfn_ndgrid import MFnNDGrid, GridPoint
from math import pi, sqrt
from numpy import sin, cos, c_, arange, hstack, array, loadtxt

from time import time
from os.path import join

#-------------------------
# calculation parameters:
#-------------------------

#--- discretization ---:
# 'shape[0]' must be multiple of two for symmetric loading
# or a multiple of four in case of asymmetric loading
shape = (4 * 20, 2)

#--- loading type   ---:
#loading_type = 'symmetric'
loading_type = 'asymmetric'

#--- element type   ---:
#element_type = 'linear'
element_type = 'quadratic'

#--- TLine spec     ---:
# Define the value of the load and the step size:
# (Note: in the boundary conditions imposed displacements of value 1 are defined)
u_max = -0.020 # maximum imposed displacement in [m]

t_max = 1.0
n_steps = 10
tline = TLine(min=0.0, step=t_max / n_steps, max=t_max)

#time_function = lambda t: sqrt(t)

#-------------------------
# geometry:
#-------------------------

# parameters used in the geometry transformation
# parameters correspond to the slices of the cylindric shell project
beta_degree = 100.18 # [degree]
D_top_ = 0.03   # [m]
D_bottom_ = 0.04   # [m]

# the thickness in out of plane direction (z-direction)
# is captured by multiplication of thickness with the
# elasticity matrix: 
thickness = 0.35   # [m]


#-------------------------
# material model:
#-------------------------

#mp = MATS2DElastic(E = 30000.,
#                   nu = 0.2,
#                   stress_state = 'plane_stress')

#mp = MATS2DScalarDamage(E = 30000.,
#                        nu = 0.2,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              


#                        epsilon_f = .5e-2,
#                        strain_norm_type = 'Mazars')

phi_fn_type = 'strain_hardening'

if phi_fn_type == 'general':
    phi_fn = PhiFnGeneral()
    # used data of phi_fn fitted from tensile test: 
    fitted_phi_fn = loadtxt(join('mats_lab', 'fitted_phi_fn.out'))
    x = fitted_phi_fn[0]
    y = fitted_phi_fn[1]
    phi_fn.mfn.set(xdata=x, ydata=y)
    phi_fn.mfn.data_changed = True
elif phi_fn_type == 'strain_hardening':
    phi_fn = PhiFnStrainHardening(Epp=5e-5,
                                   Efp=5e-4,
                                   Dfp=0.2,
                                   Elimit=0.02
                                   )

# --- cmdm ---:
mp = MATS2DMicroplaneDamage(E=34000.*thickness,
                          nu=0.2,
                          symmetrization='sum-type',
                          model_version='compliance',
                          phi_fn=phi_fn)
#@todo: the number of microplanes can not be changed;
# an index error occurres!?
###mp.n_mp = 30

#-------------------------
# element type:
#-------------------------

if element_type == 'linear':
    fets_eval = FETS2D4Q(mats_eval=mp)
elif element_type == 'quadratic':
    fets_eval = FETS2D9Q(mats_eval=mp)


#-------------------------
# Geometry transformation:
#-------------------------

def arch_2d(points):
    x = points[:, 0]
    y = points[:, 1]

    # outer radius of the arc
    R = 1.395
    # opening angle of the arc (=segment of a circle)
    beta = beta_degree * pi / 180
    #thickness at the top of the arc
    D_top = D_top_
    #thickness at the bottom of the arc
    D_bottom = D_bottom_
    # thickness change between top and bottom
    D_delta = D_bottom - D_top
    # angle between horizontal line and starting point of the arc
    alpha = (pi - beta) / 2

    # derive the length and heights of the carthesian grid
    L = x[-1] - x[0]
    H = y[-1] - y[0]

    # --------------------------------------------------
    # transformation from carthesian coordinates (x,y) 
    # to curvilinear (polar) coordinates (phi,r):
    # --------------------------------------------------

    # angle starts at the horizontal line turning right:
    phi = alpha + x / L * beta

    # variable thickness:

#    # (linear transition assumed)
#    # if the value of phi <= (pi/2):
#    D_var = D_bottom - (phi - alpha)/(beta/2.) * D_delta    
#    # if the value of phi > (pi/2):
#    bool_phi = phi > pi/2
#    D_var[bool_phi] = D_top + (phi[bool_phi] - pi/2.)/(beta/2.) * D_delta    

    # (quadratic transition assumed)
    D_var = D_delta * ((phi - alpha) / (pi / 2 - alpha)) ** 2 - 2.0 * D_delta * ((phi - alpha) / (pi / 2 - alpha)) + D_bottom

    # radius: 
    r = R - (H - y) / H * D_var

    # carthesian coordinates transformed in the arc geometry
    x, y = -r * cos(phi), r * sin(phi)
    return c_[ x, y ]


#-------------------------
# Discretization:
#-------------------------

# The 'coord_max'-coordinates of FEGrid
# are derived in the method arch2d and 
# automatically considered in the transformation
# therefore the values define here do not effect 
# the calculation
length = 1.
height = 1.

domain = FEGrid(coord_max=(length, height, 0.),
                 shape=shape,
                 fets_eval=fets_eval,
                 geo_transform=arch_2d)


#-------------------------
# ts: 
#-------------------------

# Depends on the choice of loading specified in the calculation 
# parameters defined above ('symmetric' or 'asymmetric' loading):

# --- (single displacement at top of the arc)    
if loading_type == 'symmetric':
    # get dofs corresponding to the loading case:
    top_middle_dofs, top_middle_dof_points = domain.get_top_middle_dofs()
    # dof in global y-direction:
    top_middle_dof = top_middle_dofs[0, 1]

    # boundary condition for symmetric loading case:
    BCDOF_list = [ BCDofGroup(var='u', value= -1., dims=[1],
                          get_dof_method=domain.get_top_middle_dofs)]

    # rtrace for symmetric loading case:
    var_y_ = 'F_int'
    idx_y_ = top_middle_dofs[0, 1]
    var_x_ = 'U_k'
    idx_x_ = top_middle_dofs[0, 1]

# --- (single displacement at the left side at L/4)    
elif loading_type == 'asymmetric':
    # get dofs corresponding to the loading case:
    top_dofs, top_dof_points = domain.get_top_dofs()

    # dof in global y-direction:
    if fets_eval.n_e_dofs == 8:
        top_left_middle_dofs = top_dofs[shape[0] / 4, :]
    elif fets_eval.n_e_dofs == 18:
        top_left_middle_dofs = top_dofs[shape[0] / 2, :]

    # loading angle for radial (asymmetric) loading:
    beta = beta_degree * pi / 180
    alpha = (pi - beta) / 2
    phi_left_middle = alpha + beta / 4

    # boundary condition for asymmetric loading case
#    BCDOF_x = BCDof( var='u', value =  cos(phi_left_middle), dof = top_left_middle_dofs[0] )
#    BCDOF_y = BCDof( var='u', value = -sin(phi_left_middle), dof = top_left_middle_dofs[1] )
#    BCDOF_list = [BCDOF_x, BCDOF_y]
    BCDOF_radial = BCDof(var='u',
                          value=cos(phi_left_middle), dof=top_left_middle_dofs[0],
                          link_dofs=[top_left_middle_dofs[1]],
                          link_coeffs=[ sin(phi_left_middle) / cos(phi_left_middle)])
    BCDOF_y = BCDof(var='u', value=u_max, dof=top_left_middle_dofs[1],
                     time_function=lambda t: sqrt(t))
    BCDOF_list = [BCDOF_y]

    # rtrace for asymmetric loading case:
    var_y_ = 'F_int'
    idx_y_ = top_left_middle_dofs[1]
    var_x_ = 'U_k'
    idx_x_ = top_left_middle_dofs[1]
    U_x_ = 'U_k'
    idx_U_x_ = top_left_middle_dofs[0]
    U_y_ = 'U_k'
    idx_U_y_ = top_left_middle_dofs[1]

bcond_list = [ # constraints at supports 

#                  # in case of a clamped support:
#                  BCDofGroup( var='u', value = 0., dims = [0,1],
#                              get_dof_method = domain.get_left_dofs ),
#                  BCDofGroup( var='u', value = 0., dims = [0,1],
#                              get_dof_method = domain.get_right_dofs ),

                  # in case of a hinge support (free rotations):
                  BCDofGroup(var='u', value=0., dims=[0, 1],
                              get_dof_method=domain.get_bottom_left_dofs),
                  BCDofGroup(var='u', value=0., dims=[0, 1],
                              get_dof_method=domain.get_bottom_right_dofs),


                # 'BCDOF_list' is defined above based on specified parameter 'loading_type'
               ] + BCDOF_list

tstepper = TS(
     sdomain=domain,
     bcond_list=bcond_list,
     rtrace_list=[
                 RTraceGraph(name='Fi over u at the top of the arc (iteration)' ,
                             # defined above based on parameter 'loading_type'
                             var_x=var_x_,
                             var_y=var_y_,
                             idx_x=idx_x_,
                             idx_y=idx_y_,
                             # define a function that transforms 
                             # the rtace data for the visulization
                             transform_x='-x',
                             transform_y='-y',
                             #
                             record_on='update'),

                 RTraceGraph(name='u_y over u_x at the top of the arc (iteration)' ,
                             # defined above based on parameter 'loading_type'
                             var_x=U_x_,
                             var_y=U_y_,
                             idx_x=idx_U_x_,
                             idx_y=idx_U_y_,
                             #
                             record_on='update'),

              RTraceDomainListField(name='Displacement' ,
              var='u', idx=0,
              record_on='update',
              warp=True),

              RTraceDomainListField(name='Strain' ,
              var='eps_app', idx=0,
              record_on='update',
              warp=True),

              RTraceDomainListField(name='Stress' ,
              var='sig_app', idx=0,
              record_on='update',
              warp=True),

              # only valid for scalar damage model
              RTraceDomainListField(name='Fracture Energy' ,
              var='fracture_energy', idx=0,
              record_on='update',
              warp=True),
            ]
        )


# Add the time-loop control
tloop = TLoop(tstepper=tstepper,
               tolerance=1e-4,
               tline=tline)

tloop.eval()

# Put the whole thing into the simulation-framework to map the
# individual pieces of definition into the user interface.
#
from ibvpy.plugins.ibvpy_app import IBVPyApp
app = IBVPyApp(ibv_resource=tloop)
app.main()
