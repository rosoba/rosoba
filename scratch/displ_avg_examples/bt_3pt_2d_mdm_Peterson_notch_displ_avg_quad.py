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

from ibvpy.cntl.displ_avg.rt_nonlocal_averaging import \
    RTNonlocalAvg, QuarticAF

def app():
    avg_radius = 95
    thickness = 50.0  # mm
    md = MATS2DScalarDamage(E=30.0e3,
                            nu=0.2,
                            epsilon_0=1.0e-4,
                            epsilon_f=8.0e-4,
                            # epsilon_f = 12.0e-4, #test doubling the e_f
                            stress_state="plane_stress",
                            stiffness="secant",
                            # stiffness  = "algorithmic",
                            strain_norm=Rankine())

    mdm = MATS2DMicroplaneDamage(E=30.0e3,
                                 nu=0.2,
                                 # elastic_debug=True,
                                model_version='compliance',
                                phi_fn=PhiFnStrainSoftening(
                                                              G_f=0.124,
                                                              f_t=3.3,
                                                              md=0.0,
                                                              h=12.5))

    # mdm = MATS2DElastic()

    fets_eval = FETS2D4Q9U(mats_eval=mdm)

    fe_domain = FEDomain()

    fe_rgrid = FERefinementGrid(name='fe_grid1', fets_eval=fets_eval, domain=fe_domain)

    n_half = 10
    n_el_x = n_half * 2 + 1
    n_el_y = n_el_x / 10
    # Discretization
    fe_grid = FEGrid(coord_max=(2000.0, 200.0, 0.),
                      shape=(n_el_x, n_el_y),
                      fets_eval=fets_eval,
                      level=fe_rgrid)

    for i in range(0, n_el_y / 2):
        fe_grid.deactivate((n_half, i))

    mf = MFnLineArray(xdata=array([0, 1]),
                       ydata=array([0, 1]))

    # averaging function
#     avg_processor = RTNonlocalAvg(avg_fn=QuarticAF(radius=avg_radius,
#                                                         correction=True))

    loading_dof = fe_grid[n_el_x / 2, -1, 0, -1].dofs.flatten()[1]
    loading_dofs = fe_grid[n_el_x / 2, -1, (0, -1), -1].dofs.flatten()
    print 'loading_dofs', loading_dofs

    redundant_dofs = fe_grid[n_el_x / 2., :n_el_y / 2, 1, :-1]
    print 'redundant_dof' , redundant_dofs.dofs.flatten()

    ts = TS(sdomain=fe_grid,
            # u_processor=avg_processor,
             bcond_list=[
                        # constraint for all left dofs in y-direction:
                        BCSlice(var='u', slice=fe_grid[0, 0, 0, 0], dims=[0, 1], value=0.),
                        BCSlice(var='u', slice=fe_grid[-1, 0, -1, 0], dims=[1], value=0.),
                        BCSlice(var='u', slice=redundant_dofs, dims=[0, 1], value=0.),
                        BCSlice(var='u', slice=fe_grid[n_el_x / 2, -1, (0, -1), -1], dims=[1],
                                time_function=mf.get_value,
                                value= -1.0),

                        ],
             rtrace_list=[RTraceGraph(name='Fi,right over u_right (iteration)' ,
                                      var_y='F_int', idx_y_arr=loading_dofs,
                                      var_x='U_k', idx_x=loading_dof,
                                      transform_x='-x', transform_y='-y*%g' % thickness,
                                      record_on='update'),
                          RTraceDomainListField(name='Strain' ,
                                      var='eps_app', idx=0,
                                      record_on='update'),
                          RTraceDomainListField(name='Stress' ,
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

    #                        RTraceDomainField(name = 'Stress' ,
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

    tl.eval()
#    # Put the whole stuff into the simulation-framework to map the
#    # individual pieces of definition into the user interface.
#    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    ibvpy_app = IBVPyApp(ibv_resource=ts)
    ibvpy_app.main()


if __name__ == '__main__':
    app()
