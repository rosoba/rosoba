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
    avg_radius = 0.1

    md = MATS2DScalarDamage(E=20.0e3,
                            nu=0.2,
                            epsilon_0=1.0e-4,
                            epsilon_f=8.0e-4,
                            # epsilon_f = 12.0e-4, #test doubling the e_f
                            stress_state="plane_stress",
                            stiffness="secant",
                            # stiffness  = "algorithmic",
                            strain_norm=Rankine())

    mdm = MATS2DMicroplaneDamage(E=20.0e3,
                                 nu=0.2,
                                 # epsilon_f = 12.0e-4, #test doubling the e_f
                    #    tl.eval()
#    # Put the whole stuff into the simulation-framework to map the
#    # individual pieces of definition into the user interface.
#    #
#    from ibvpy.plugins.ibvpy_app import IBVPyApp
#    ibvpy_app = IBVPyApp(ibv_resource=ts)
#    ibvpy_app.main()             stress_state="plane_strain",
                                model_version='compliance',
                                phi_fn=PhiFnStrainSoftening(
                                                              G_f=0.0014,
                                                              f_t=2.0,
                                                              md=0.0,
                                                              h=2. * avg_radius))
    print 'Microplane normals' , mdm._get__MPN()
    print 'Microplane weights' , mdm._get__MPW()

    fets_eval = FETS2D4Q(mats_eval=mdm)  # , ngp_r = 3, ngp_s = 3)

    n_el_x = 20
    # Discretization
    fe_grid = FEGrid(coord_max=(2, .2, 0.),
                      shape=(n_el_x, n_el_x / 4),
                      fets_eval=fets_eval)

    mf = MFnLineArray(xdata=array([0, 1, 2 ]),
                       ydata=array([0, 3., 3.2 ]))

    # averaging function
    avg_processor = RTNonlocalAvg(avg_fn=QuarticAF(radius=avg_radius,
                                                       correction=True))

    loading_dof = fe_grid[n_el_x / 2, -1, 0, -1].dofs.flatten()[1]
    print 'loading_dof', loading_dof
    ts = TS(sdomain=fe_grid,
             u_processor=avg_processor,
             bcond_list=[
                        # constraint for all left dofs in y-direction:
                        BCSlice(var='u', slice=fe_grid[0, 0, 0, 0], dims=[0, 1], value=0.),
                        BCSlice(var='u', slice=fe_grid[-1, 0, -1, 0], dims=[1], value=0.),
                        BCSlice(var='u', slice=fe_grid[n_el_x / 2, -1, 0, -1], dims=[1],
                                time_function=mf.get_value,
                                value= -2.0e-4),
                        ],
             rtrace_list=[
                            RTraceGraph(name='Fi,right over u_right (iteration)' ,
                                      var_y='F_int', idx_y=loading_dof,
                                      var_x='U_k', idx_x=loading_dof,
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
                tline=TLine(min=0.0, step=.1, max=1.0))

    tl.setup()
    print avg_processor.C_mtx


    tl.eval()
#    # Put the whole stuff into the simulation-framework to map the
#    # individual pieces of definition into the user interface.
#    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    ibvpy_app = IBVPyApp(ibv_resource=ts)
    ibvpy_app.main()


if __name__ == '__main__':
    app()
