"""Fabric (fabfile.org) tasks for SMC."""

from fabric.api import *
from pyboxlib.fabutils import *

@task
def flameball():
    """Convergence tests for the FlameBall example."""

    env.probin = 'inputs-cflconv'
    env.exe    = 'main.*.exe'

    convergence('cflconv', 'cflfac', [ 0.1, 0.2, 0.4, 0.6, 0.8, 1.0 ], 
                sdc_nnodes=5, sdc_iters=10)
