'''
Created on Jul 21, 2010

@author: jakub
'''
from ibvpy.api import \
    TStepper as TS, RTraceGraph, TLoop, \
    TLine, BCDof, IBVPSolve as IS, DOTSEval, FEDomain, FERefinementGrid, \
    FEGrid, BCSlice, RTraceDomainListField


# from apps.scratch.jakub.mlab.mlab_trace import RTraceDomainListField
from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import MATS2DScalarDamage
from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm import MATS2DMicroplaneDamage, PhiFnStrainSoftening
from ibvpy.mats.mats2D.mats2D_sdamage.strain_norm2d import Euclidean, Mazars, Rankine
from ibvpy.mats.mats_proxy import MATSProxy
from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
from ibvpy.fets.fets2D.fets2D4q9u import  FETS2D4Q9U
from ibvpy.fets.fets2D.fets2D4q8u import FETS2D4Q8U
from numpy import array, cos, sin, pi, sqrt, deg2rad, arctan
from mathkit.mfn import MFnLineArray

from rt_nonlocal_averaging import \
    RTNonlocalAvg, QuarticAF

def app():
    avg_radius = 30

    md = MATS2DScalarDamage(E=30.0e3,
                            nu=0.2,
                            epsilon_0=1.0e-4,
                            epsilon_f=8.0e-4,
                            # epsilon_f = 12.0e-4, #test doubling the e_f
                            stress_state="plane_strain",
                            stiffness="secant",
                            # stiffness  = "algorithmic",
                            strain_norm=Rankine())

    mdm = MATS2DMicroplaneDamage(E=30.0e3,
                                 nu=0.2,
                                 # epsilon_f = 12.0e-4, #test doubling the e_f
                                 stress_state="plane_stress",
                                 model_version='compliance',
                                 phi_fn=PhiFnStrainSoftening(
                                                              G_f=0.14,
                                                              f_t=3.3,
                                                              md=0.0,
                                                              h=2. * avg_radius))

    fets_eval = FETS2D4Q(mats_eval=mdm)  # , ngp_r = 3, ngp_s = 3)

    n_el_x = 2
    # Discretization
    fe_grid = FEGrid(coord_max=(40, 20, 0.),
                      shape=(n_el_x, n_el_x / 2),
                      fets_eval=fets_eval)

    mf = MFnLineArray(xdata=array([0, 1]),
                       ydata=array([0, 1]))

    Ps_fn = MFnLineArray(xdata=array([0, 3, 10 ]),
                       ydata=array([0, 1, 1 ]))

    P_fn = MFnLineArray(xdata=array([0, 3, 10 ]),
                       ydata=array([0, 0, 1 ]))

    # averaging function
#     avg_processor = RTNonlocalAvg(avg_fn=QuarticAF(radius=avg_radius,
#                                                        correction=True))


    loading_slice = fe_grid[-1 , -1, -1, -1]  # left top Node 3 [2,3]
#     loading_slice = fe_grid[n_el_x - 1 , -1, -1, -1] #right top Node 2 [6,7]
#     loading_slice = fe_grid[n_el_x - 1 , -1, -1, 0]  # right down Node 1 [4,5]
#     loading_slice = fe_grid[n_el_x - 1 , -1, 0, 0] #left down Node 0 [0,1]


    loading_dof = loading_slice.dofs.flatten()[1]
    loading_dofs = loading_slice.dofs.flatten()
    print 'loading_dofs', loading_dofs

    print 'loading_dof', loading_dof
    support_dof_left = fe_grid[0, 0, 0, 0].dofs[0, 0, 1]
    support_dof_right = fe_grid[-1, 0, -1, 0].dofs[0, 0, 1]
    print 'support_dof', support_dof_left, support_dof_right

    Linked_Slice = fe_grid[0, 0, -1, -1].dofs.flatten()[1]
    Linked_Dofs = (fe_grid[0, 0, 0, -1].dofs.flatten()[1], fe_grid[-1, 0, -1, -1].dofs.flatten()[1])

    print 'Linked_Slice', Linked_Slice
    print 'Linked_Dofs', Linked_Dofs

    dof_ref = fe_grid[-1, -1, -1, -1].dofs[0, 0, 1].flatten()[0]
    print 'dof_ref', dof_ref



    link_up = BCSlice(var='u', value=0., dims=[1],
                                   slice=fe_grid[:, 0, :, -1],
                                   link_slice=fe_grid[-1, -1, -1, -1],
                                   link_dims=[1],
                                   link_coeffs=[1.])

    ts = TS(sdomain=fe_grid,
#              u_processor=avg_processor,
             bcond_list=[
                        link_up,
                        # constraint for all left dofs in y-direction:
#                         BCSlice(var='u', slice=fe_grid[0, 0, 0, -1], link_slice=fe_grid[0, 0, -1, -1], link_coeffs=[-1, 1]),
                        BCSlice(var='u', slice=fe_grid[0, 0, 0, 0], dims=[0, 1], value=0.),
                        BCSlice(var='u', slice=fe_grid[-1, 0, -1, 0], dims=[1], value=0.),
                        BCSlice(var='f', slice=loading_slice, dims=[1],
                                time_function=Ps_fn.get_value,
                                value= -50),
                        ],
             rtrace_list=[
                            RTraceGraph(name='F(loading point) - w' ,
                                        var_y='F_int', idx_y_arr=loading_dofs,
                                        transform_y='-y',
                                        var_x='U_k', idx_x=loading_dof,
                                        transform_x='-x',
                                        record_on='update'),
                            RTraceGraph(name='F(support point) - w' ,
                                      var_y='F_int', idx_y_arr=[support_dof_left, support_dof_right],
                                      transform_y='y',
                                      var_x='U_k', idx_x=loading_dof,
                                      transform_x='-x',
                                      record_on='update'),
                            RTraceDomainListField(name='Deformation' ,
                                           var='eps_app', idx=0,
                                           record_on='update'),
                            RTraceDomainListField(name='Deformation' ,
                                           var='sig_app', idx=0,
                                           record_on='update'),
                            RTraceDomainListField(name='Displacement' ,
                                           var='u', idx=1,
                                           record_on='update',
                                           warp=True),
                            RTraceDomainListField(name='fracture_energy' ,
                                           var='fracture_energy', idx=0,
                                           record_on='update',
                                           warp=True),
                            RTraceDomainListField(name='Damage' ,
                                        var='omega_mtx', idx=0, warp=True,
                                        record_on='update'),

    #                         RTraceDomainField(name = 'Stress' ,
    #                                        var = 'sig', idx = 0,
    #                                        record_on = 'update'),
    #                        RTraceDomainField(name = 'N0' ,
    #                                       var = 'N_mtx', idx = 0,
    #                                       record_on = 'update')
                        ]
            )

    # Add the time-loop control
    #
    tl = TLoop(tstepper=ts,
                tolerance=5.0e-5,
                KMAX=200,
                tline=TLine(min=0.0, step=.5, max=10))

    tl.setup()
#     print avg_processor.C_mtx


    tl.eval()
#    # Put the whole stuff into the simulation-framework to map the
#    # individual pieces of definition into the user interface.
#    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    ibvpy_app = IBVPyApp(ibv_resource=ts)
    ibvpy_app.main()


if __name__ == '__main__':
    app()
