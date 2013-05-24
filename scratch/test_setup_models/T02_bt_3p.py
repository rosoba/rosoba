'''
Three point bending test
'''

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
    TLine, BCDof, IBVPSolve as IS, BCDofGroup

from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
from ibvpy.fets.fets2D.fets2D9q import FETS2D9Q

from ibvpy.mesh.fe_grid import FEGrid

from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import \
    MATS2DElastic
from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import \
    MATS2DScalarDamage
from ibvpy.mats.mats2D.mats2D_sdamage.strain_norm2d import \
    Energy, Euclidean, Mises, Rankine, Mazars
from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm import \
    MATS2DMicroplaneDamage, \
    PhiFnGeneral, PhiFnStrainSoftening, PhiFnStrainHardening
from ibvpy.mats.mats_proxy import \
    MATSProxy

from mathkit.geo.geo_ndgrid import \
    GeoNDGrid
from mathkit.mfn.mfn_ndgrid.mfn_ndgrid import \
    MFnNDGrid, GridPoint

from math import \
    pi as Pi, cos, sin, exp, sqrt as scalar_sqrt

from numpy import \
    sin, cos, c_, arange, hstack, array, loadtxt

from time import time

from os.path import join


#-------------------------
# calculation parameters:
#-------------------------

# geometry:
length = 0.80
height = 0.04
thickness = 0.15

# discretization:
# Note: For loading at the center at the top an even 
#       number of elements in x-direction is necessary.
shape = (30, 4)

# loading:
u_max = -0.0017
t_max = 1.0
n_steps = 16
tline = TLine(min=0.0, step=t_max / n_steps, max=t_max)

#element_type = 'linear'
element_type = 'quadratic'

reinforcement = True
isotropy = False


#-------------------------
# material model:
#-------------------------

cmdm = MATS2DMicroplaneDamage(E=34000. * thickness,
                              nu=0.2,
                              symmetrization='sum-type',
                              model_version='compliance',
                              )
#--- unreinforced ---
if reinforcement == False:
    # Isotropic, QuasiBrittle:
    cmdm.phi_fn = PhiFnStrainSoftening()

#--- reinforced ---
elif reinforcement == True and isotropy == False:
#    ## Anisotropic, TensionStiffening:
    cmdm.phi_fn = PhiFnStrainHardening(Epp=59e-6,
                                        Efp=191e-6,
                                        Dfp=0.4,
                                        Elimit=1.e-3)
    # 'polar_fn': parameter definition: 
    cmdm.varied_params = ['Dfp']
    cmdm.varpars['Dfp'].polar_fn.set(phi_residual=1., \
                                     phi_quasibrittle=0.0,
                                     alpha=0,
                                     delta_trans=Pi / 4. ,
                                     delta_alpha=Pi / 8.)

elif reinforcement == True and isotropy == True:

    # used data of phi_fn fitted from tensile test: 
    fitted_phi_fn = loadtxt(join('mats_lab', 'fitted_phi_fn.out'))
    x = fitted_phi_fn[0]
    y = fitted_phi_fn[1]

    cmdm.phi_fn = PhiFnGeneral()
    cmdm.phi_fn.mfn.set(xdata=x, ydata=y)
    cmdm.phi_fn.mfn.data_changed = True

#-------------------------
# element type:
#-------------------------
if element_type == 'linear':
    fets_eval = FETS2D4Q(mats_eval=cmdm)
elif element_type == 'quadratic':
    fets_eval = FETS2D9Q(mats_eval=cmdm)


#-------------------------
# Discretization
#-------------------------

domain = FEGrid(coord_max=(length, height, 0.),
                 shape=shape,
                 fets_eval=fets_eval)


#-------------------------
# ts: 
#-------------------------

# alternative way to get the top/bottom middle dof using slicing of FEGrid:
# i.e. get the element which is placed right for the center and then get the lower left node

#bottom_middle_dof = domain[shape[0]/2, 0, 0, 0].dofs[0,0,1]
#top_middle_dof    = domain[shape[0]/2,-1, 0,-1].dofs[0,0,1]

top_middle_dofs, top_middle_dof_points = domain.get_top_middle_dofs()
top_middle_dof = top_middle_dofs[0, 1]
print 'top_middle_dof', top_middle_dof

bottom_middle_dofs, bottom_middle_dof_points = domain.get_bottom_middle_dofs()
bottom_middle_dof = bottom_middle_dofs[0, 1]
print 'bottom_middle_dof', bottom_middle_dof

tstepper = TS(dof_resultants=True,
               sdomain=domain,
     bcond_list=[ BCDof(var='u', value=u_max, dof=top_middle_dof),
                     BCDofGroup(var='u', value=0., dims=[0, 1],
                              get_dof_method=domain.get_bottom_left_dofs),
                     BCDofGroup(var='u', value=0., dims=[1],
                              get_dof_method=domain.get_bottom_right_dofs) ],
     rtrace_list=[
                 RTraceGraph(name='Fi at top over u at bottom (iteration)' ,
                             var_y='F_int', idx_y=top_middle_dof,
                             var_x='U_k', idx_x=bottom_middle_dof,
                             transform_x='-x',
                             transform_y='-y',
                             record_on='update'
                             ),
              RTraceDomainListField(name='Displacement' ,
              var='u', idx=0,
              record_on='update',
              warp=True),

              RTraceDomainListField(name='Strain' ,
              var='eps_app', idx=0,
              record_on='update',
              warp=False),

              RTraceDomainListField(name='Stress' ,
              var='sig_app', idx=0,
              record_on='update',
              warp=False),

              RTraceDomainListField(name='Fracture Energy' ,
              var='fracture_energy', idx=0,
              record_on='update',
              warp=True),
            ]
        )

# Add the time-loop control
tloop = TLoop(tstepper=tstepper, tolerance=1e-4,
               tline=tline)

tloop.eval()

# Put the whole thing into the simulation-framework to map the
# individual pieces of definition into the user interface.
#
from ibvpy.plugins.ibvpy_app import IBVPyApp
app = IBVPyApp(ibv_resource=tloop)
app.main()

