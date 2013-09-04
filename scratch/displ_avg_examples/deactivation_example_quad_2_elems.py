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

    mdm = MATS2DElastic()

    fets_eval = FETS2D4Q(mats_eval=mdm)

    fe_domain = FEDomain()

    fe_rgrid = FERefinementGrid(name='fe_grid1', fets_eval=fets_eval, domain=fe_domain)

    fe_grid = FEGrid(coord_max=(2.0, 1.0),
                      shape=(2, 1),
                      fets_eval=fets_eval,
                      level=fe_rgrid)

    for i in range(0, 1):
        fe_grid.deactivate((0, 0))

    ts = TS(sdomain=fe_grid,
             bcond_list=[
                         BCSlice(var='u', slice=fe_grid[-1, :, -1, :], dims=[0, 1], value=0.),
                         BCSlice(var='u', slice=fe_grid[0, 0, :-1, :], dims=[0, 1], value=0.),
                         BCSlice(var='u', slice=fe_grid[-1, 0, 0, -1], dims=[1],
                                 value= -1.0),
                        ],
            rtrace_list=[
                          RTraceDomainListField(name='Strain' ,
                                      var='eps_app', idx=0,
                                      record_on='update'),
                          RTraceDomainListField(name='Stress' ,
                                      var='sig_app', idx=0,
                                      record_on='update'),
                          RTraceDomainListField(name='Displacement' ,
                                      var='u', idx=1,
                                      record_on='update',
                                      warp=True)
                         ]
            )

    tl = TLoop(tstepper=ts,
                tolerance=5.0e-5,
                KMAX=200,
                tline=TLine(min=0.0, step=.1, max=0.1))

    tl.eval()

    from ibvpy.plugins.ibvpy_app import IBVPyApp
    ibvpy_app = IBVPyApp(ibv_resource=ts)
    ibvpy_app.main()


if __name__ == '__main__':
    app()
