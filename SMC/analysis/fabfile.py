"""Fabric (fabfile.org) tasks for SMC."""

from math import ceil

from fabric.api import *
from pyboxlib.utils import *
from itertools import product

from fabsubmit import *
from fabutils import *

projects = '/home/memmett/projects/'
work     = '/scratch/scratchdirs/memmett/'

env.rsync = [ (projects + 'Combustion/SMC', work + 'Combustion/SMC') ]
env.bin   = work + 'Combustion/SMC/bin/FlameBall'


def setenv():
  if env.host == 'hopper':
    env.host_string = 'hopper.nersc.gov'
    env.exe = 'main.Linux.Cray.mpi.omp.exe'

def walltime(seconds):
  import math

  s = int(math.ceil(seconds))
  h = s / 3600
  m = (s % 3600) / 60
  s = s - h*3600 - m*60

  return "%02d:%02d:%02d" % (h, m, s)


@task
def flameball_cflconv():
  """Convergence tests for the FlameBall example.

  This task should be run in the bin/FlameBall directory.
  """
  setenv()

  env.nthreads = 6
  env.nprocs   = 4
  env.pernode  = 4

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

      runs.append(('sdcgl%d_%.2f' % (nnodes, cflfac), rundir, walltime(nnodes, cflfac)))


  for name, rundir, wall in runs:
    submit(name, rundir, 'probin.nml', walltime=wall)


@task
def flameball_stconv():
  """Convergence tests for the FlameBall example."""
  setenv()

  feval = 0.01             # seconds per function eval on 32**3 domain

  env.nthreads = 6
  env.nprocs   = 6
  env.pernode  = 4

  probin = Probin('inputs-stconv')
  runs   = []

  with stage() as s:

    nx0 = 32
    dt0 = 1.0e-7

    for nx, cfl, nnodes in product( [ 32, 64, 128 ],
                                    [ 0.25, 0.5, 0.75, 1.0 ],
                                    [ 3, 5 ] ):

      stop_time = 1.0e-5
      dt        = dt0 * float(nx0) / nx * cfl
      nsteps    = int(ceil(stop_time / dt))

      max_grid_size = 32
      nprocs        = nx**3 / max_grid_size**3

      name   = 'nx%03d_gl%d_cfl%.2f' % (nx, nnodes, cfl)
      rundir = 'stconv/' + name
      spath  = s.mkdir(rundir)

      probin.update(
          nx=nx, dt=dt, stop_time=stop_time,
          sdc_nnodes=nnodes,
          sdc_iters=2*nnodes-2,
          max_grid_size=max_grid_size,
          advance_method=2)
      probin.write(spath + 'probin.nml')

      wtime = walltime(max(feval * (nx/32)**3 * nnodes * (2*nnodes-2) * nsteps, 600))

      runs.append(dict(
          name=name, rundir=rundir, nprocs=nprocs,
          probin='probin.nml', walltime=wtime))

  for run in runs:
    submit(**run)
