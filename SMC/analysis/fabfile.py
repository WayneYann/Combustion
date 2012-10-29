"""Fabric (fabfile.org) tasks for SMC."""

from fabric.api import *
from pyboxlib.fabutils import *
from pyboxlib.utils import *
from fabsubmit import *
from itertools import product


env.rsync_dirs = [ '~/projects/BoxLib/Src/Python', '~/projects/Combustion/SMC' ]
env.bin = '~/projects/Combustion/SMC/bin/FlameBall'


def setenv():
  if env.host == 'hopper':
    env.host_string = 'hopper.nersc.gov'
    env.exe = 'main.Linux.Cray.mpi.omp.exe'


@task
def flameball():
    """Convergence tests for the FlameBall example.

    This task should be run in the bin/FlameBall directory.
    """
    setenv()

    env.nthreads = 4
    env.nprocs   = 6
    env.pernode  = 3

    #rsync()
    #make()

    probin = Probin('inputs-cflconv')
    runs   = []

    with stage() as s:

        for cflfac in [ 0.1, 0.2, 0.4, 0.6 ]:
            rundir = 'cflconv/rk3_cflfac%f' % (cflfac,)
            spath  = s.mkdir(rundir)

            probin.update(sdc_nnodes=1,
                          sdc_iters=1,
                          cflfac=cflfac,
                          advance_method=1)
            probin.write(spath + 'probin.nml')

            runs.append(('rk3%f' % (cflfac,), rundir))


        for nnodes, cflfac in product([ 3, 5, 9 ], [ 0.1, 0.2, 0.4, 0.6, 0.8, 1.0 ]):
            rundir = 'cflconv/nnodes%d_cflfac%f' % (nnodes, cflfac)
            spath  = s.mkdir(rundir)

            probin.update(sdc_nnodes=nnodes,
                          sdc_iters=2*nnodes-1,
                          cflfac=cflfac,
                          advance_method=2)
            probin.write(spath + 'probin.nml')

            runs.append(('sdc%d%f' % (nnodes, cflfac), rundir))

    for name, rundir in runs:
        submit(name, rundir, 'probin.nml')
