"""Fabric (fabfile.org) tasks for SMC."""

import numpy as np
import stconv

from fabric.api import *
from jobtools import JobQueue, Job


@task
def flameball_mrconv():
  """Convergence tests for the FlameBall example using MRSDC."""

  setenv()
  #make()

  jobs = JobQueue(rwd=env.rwd, queue='regular')
  dt0  = 5e-9

  nx = 32

  # reference run

  dt     = dt0/8
  nnodes = 3
  nprocs = 1

  name = 'nx%03d_gl%d_dt%g' % (nx, nnodes, dt)
  job  = Job(name=name, param_file='inputs-mrconv', rwd='mrconv/'+name, width=nprocs)

  job.update_params(
    nx=nx,
    fixed_dt=dt,
    stop_time=dt0*10,
    cflfac=-1.0,
    sdc_nnodes=nnodes,
    sdc_nnodes_chemistry=-1,
    sdc_iters=2*nnodes-2,
    max_grid_size=nx,
    advance_method=2)

  # jobs.add(job)


  # multi-rate runs

  for dt in [ 1.5*dt0, dt0/2, dt0/4 ]:

    nnodes = (3, 9)
    nprocs = 1

    name = 'nx%03d_gl%d-%d_dt%g' % (nx, nnodes[0], nnodes[1], dt)
    job  = Job(name=name, param_file='inputs-mrconv', rwd='mrconv/'+name, width=nprocs)

    job.update_params(
      nx=nx,
      fixed_dt=dt,
      stop_time=dt0*10,
      cflfac=-1.0,
      sdc_nnodes=nnodes[0],
      sdc_nnodes_chemistry=nnodes[1],
      sdc_iters=2*nnodes[0]-2,
      max_grid_size=nx,
      advance_method=3)

    jobs.add(job)


  # single-rate runs

  for dt in [ dt0/4, dt0/8, dt0/16 ]:

    nnodes = 3
    nprocs = 1

    name = 'nx%03d_gl%d_dt%g' % (nx, nnodes, dt)
    job  = Job(name=name, param_file='inputs-mrconv', rwd='mrconv/'+name, width=nprocs)

    job.update_params(
      nx=nx,
      fixed_dt=dt,
      stop_time=dt0*10,
      cflfac=-1.0,
      sdc_nnodes=nnodes,
      sdc_nnodes_chemistry=-1,
      sdc_iters=2*nnodes-2,
      max_grid_size=nx,
      advance_method=2)

    jobs.add(job)


  # runge-kutta runs

  for dt in [ dt0/4, dt0/8, dt0/16 ]:

    nprocs = 1

    name = 'nx%03d_rk3_dt%g' % (nx, dt)
    job  = Job(name=name, param_file='inputs-mrconv', rwd='mrconv/'+name, width=nprocs)

    job.update_params(
      nx=nx,
      fixed_dt=dt,
      stop_time=dt0*10,
      cflfac=-1.0,
      sdc_nnodes=-1,
      sdc_nnodes_chemistry=-1,
      sdc_iters=-1,
      max_grid_size=nx,
      advance_method=1)

    jobs.add(job)


  jobs.submit_all()
    


@task
def flameball_stconv():
  """Convergence tests for the FlameBall example."""

  setenv()
  make()

  jobs = JobQueue(rwd=env.rwd, queue='regular')

  # convergence runs
  for nx, dt, nnodes in stconv.runs:

    max_grid_size = min(nx, 32)
    nprocs        = nx**3 / max_grid_size**3
    # max_grid_size = nx
    # nprocs        = 1

    name = 'nx%03d_gl%d_dt%g' % (nx, nnodes, dt)
    job  = Job(name=name, param_file='inputs-stconv', rwd='stconv/'+name, width=nprocs)

    wtime = '00:30:00'
    if nx == 64:
      wtime = '01:00:00'
    if nx == 128:
      wtime = '04:00:00'
    if nx == 256:
      wtime = '06:00:00'


    job.attrs['walltime'] = wtime

    job.update_params(
      nx=nx,
      fixed_dt=dt,
      stop_time=stconv.stop_time,
      cflfac=None,
      sdc_nnodes=nnodes,
      sdc_iters=2*nnodes-1,
      max_grid_size=max_grid_size,
      advance_method=2)

    jobs.add(job)

  jobs.submit_all()


@task
def flameball_stconv_comp():
  """Compute convergence errors (run stconv-comp.py) on the remote host."""
  setenv()
  rsync()

  if env.host == 'hopper':
    with prefix('module load numpy'):
      with cd(env.scratch + '/Combustion/SMC/analysis'):
        run('python stconv-comp.py')
  elif env.host[:5] == 'gigan':
    with cd(env.scratch + '/Combustion/SMC/analysis'):
      run('/home/memmett/venv/base/bin/python stconv-comp.py')


@task
def flameball_mrconv_comp():
  """Compute convergence errors (run mrconv-comp.py) on the remote host."""
  setenv()
  rsync()

  if env.host == 'hopper':
    with prefix('module load numpy'):
      with cd(env.scratch + '/Combustion/SMC/analysis'):
        run('python mrconv-comp.py')
  elif env.host[:5] == 'gigan':
    with cd(env.scratch + '/Combustion/SMC/analysis'):
      run('/home/memmett/venv/base/bin/python mrconv-comp.py')
    

@task
def flameball_stconv_plot():
  """Get convergence results and plot them (run stconv-plot.py) on the local host."""
  setenv()
  get(env.scratch + 'Combustion/SMC/analysis/stconv.pkl', 'stconv.pkl')
  local('python stconv-plot.py')


@task
def flameball_mrconv_plot():
  """Get convergence results and plot them (run mrconv-plot.py) on the local host."""
  setenv()
  get(env.scratch + 'Combustion/SMC/analysis/mrconv.pkl', 'mrconv.pkl')
  local('python mrconv-plot.py')


@task
def rsync():
  """Push (rsync) directories in env.rsync to env.host."""

  setenv()
  if env.host == 'localhost':
    return

  for src, dst in env.rsync:
      command = "rsync -avz -F {src}/ {host}:{dst}".format(
          host=env.host_rsync, src=src, dst=dst)
      local(command)


@task
def make(target=''):
  """Run make in the env.rwd directory on the remote host."""

  setenv()
  rsync()
  with cd(env.rwd):
    run('make %s' % target)


def setenv():
  """Setup Fabric and jobtools environment."""

  projects = '/home/memmett/projects/'

  if env.host == 'hopper':
    env.scratch     = '/scratch/scratchdirs/memmett/'
    env.scheduler   = 'hopper'
    env.host_string = 'hopper.nersc.gov'
    env.host_rsync  = 'hopper-s'
    env.exe         = 'main.Linux.gfortran.mpi.omp.exe'

    env.depth   = 6
    env.pernode = 4

    env.aprun_opts = [ '-cc numa_node' ]

  elif env.host[:5] == 'gigan':
    env.scratch     = '/scratch/memmett/'
    env.scheduler   = 'serial'
    env.host_string = 'gigan.lbl.gov'
    env.host_rsync  = 'gigan-s'
    env.exe         = 'main.Linux.gfortran.debug.prof.omp.exe'

    env.width = 1
    env.depth = 16
    
  env.rsync = [ (projects + 'Combustion', env.scratch + 'Combustion'),
                (projects + 'SDCLib',     env.scratch + 'SDCLib') ]

  env.rwd = env.scratch + 'Combustion/SMC/bin/FlameBall'

