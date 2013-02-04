"""Fabric (fabfile.org) tasks for SMC."""

import sys
sys.path.append('/home/memmett/projects/BoxLib/Src/Python')
sys.path.append('/home/memmett/projects/Combustion/SMC/analysis')

import numpy as np

from pyboxlib.utils import *
from fabric.api import *
from fabsubmit import *
from fabutils import *

projects = '/home/memmett/projects/'
work     = '/scratch/scratchdirs/memmett/'
#work     = '/scratch/memmett/'

env.rsync = [ (projects + 'Combustion/SMC', work + 'Combustion/SMC') ]
#env.rsync.append((projects + 'BoxLib/Src/Python', '~/projects/BoxLib/Src/Python'))
env.bin   = work + 'Combustion/SMC/bin/FlameBall'


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
  rsync()
  make()

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

      runs.append(('rk3%f' % (cflfac,), rundir, walltime(60*60*2)))

    for nnodes, cflfac in product([ 3, 5 ], [ 0.1, 0.2, 0.4, 0.6, 0.8, 1.0 ]):
      rundir = 'cflconv/nnodes%d_cflfac%f' % (nnodes, cflfac)
      spath  = s.mkdir(rundir)

      probin.update(sdc_nnodes=nnodes,
              sdc_iters=2*nnodes-1,
              cflfac=cflfac,
              advance_method=2)
      probin.write(spath + 'probin.nml')

      runs.append(('sdcgl%d_%.2f' % (nnodes, cflfac), rundir, walltime(60*60*3)))


  for name, rundir, wall in runs:
    submit(name, rundir, 'probin.nml', walltime=wall)


@task
def flameball_stconv():
  """Convergence tests for the FlameBall example."""

  import stconv

  setenv()
  # rsync()
  # make()

  env.nthreads = 6
  env.nprocs   = 6
  env.pernode  = 4

  probin = Probin('inputs-stconv')
  runs   = []

  with stage() as s:

    # convergence runs
    for nx, dt, nnodes in stconv.runs:

      max_grid_size = min(nx, 32)
      nprocs        = nx**3 / max_grid_size**3

      name   = 'nx%03d_gl%d_dt%g' % (nx, nnodes, dt)
      rundir = 'stconv/' + name
      spath  = s.mkdir(rundir)

      probin.update(
          nx=nx, fixed_dt=dt, stop_time=stconv.stop_time, cflfac=None,
          sdc_nnodes=nnodes,
          sdc_iters=2*nnodes-1,
          max_grid_size=max_grid_size,
          advance_method=2)
      probin.write(spath + 'probin.nml')

      wtime = walltime(2*60*60)

      runs.append(dict(
          name=name, rundir=rundir, nprocs=nprocs,
          probin='probin.nml', walltime=wtime))

  for run in runs:
    submit(**run)


@task
def flameball_stconv_comp():
  setenv()
  rsync()
  
  with prefix('module load numpy'):
    with cd(work + '/Combustion/SMC/analysis'):
      run('python stconv-comp.py')


@task
def flameball_stconv_plot():
  setenv()
  
  get(work + 'Combustion/SMC/analysis/stconv.pkl', 'stconv.pkl')
  local('python stconv-plot.py')
