"""Fabric (fabfile.org) tasks for SMC."""

import numpy as np
#import stconv

from fabric.api import *
from jobtools import JobQueue, Job
from itertools import product



###############################################################################
# flamebox multi-rate speed tests

@task
def flamebox_speed():
  """Timing tests for the FlameInABox example using MRSDC."""

  setenv('Combustion/SMC/bin/FlameInABox', find_exe=True)
  make()

  jobs = JobQueue(queue='regular', walltime="01:00:00")

  dt0        = 6e-9
  start_time = 0
  stop_time  = start_time + dt0 * 70

  nx            = 32
  max_grid_size = 32
  nprocs        = nx**3 / max_grid_size**3

  # reference run
  dt   = 1e-9
  name = 'gl3_dt%e' % dt

  job = Job(name=name, param_file='inputs-flamebox', rwd='speed/'+name, 
            width=nprocs, walltime="02:00:00")
  job.update_params(
    advance_method=2, stop_time=stop_time,
    sdc_nnodes=3, sdc_nnodes_chemistry=0, sdc_iters=4,
    fixed_dt=dt, nx=nx, max_grid_size=max_grid_size)
  jobs.add(job)

  # timing runs
  for dt0 in [ 6e-9, 3e-9, 1.5e-9 ]:

    # single-rate sdc run
    dt   = dt0
    name = 'gl3_dt%e' % dt

    job = Job(name=name, param_file='inputs-flamebox', rwd='speed/'+name, 
              width=nprocs, walltime="02:00:00")
    job.update_params(
      advance_method=2, stop_time=stop_time,
      sdc_nnodes=3, sdc_nnodes_chemistry=0, sdc_iters=4,
      fixed_dt=dt, nx=nx, max_grid_size=max_grid_size)
    jobs.add(job)

    # multi-rate sdc runs
    dt   = 5*dt0
    name = 'gl3.9_dt%e' % dt

    job = Job(name=name, param_file='inputs-flamebox', rwd='speed/'+name, width=nprocs)
    job.update_params(
      advance_method=3, stop_time=stop_time,
      sdc_nnodes=3, sdc_nnodes_chemistry=9, sdc_iters=5,
      fixed_dt=dt, nx=nx, max_grid_size=max_grid_size)
    jobs.add(job)

    dt   = 7*dt0
    name = 'gl3.13_dt%e' % dt

    job = Job(name=name, param_file='inputs-flamebox', rwd='speed/'+name, width=nprocs)
    job.update_params(
      advance_method=3, stop_time=stop_time,
      sdc_nnodes=3, sdc_nnodes_chemistry=13, sdc_iters=5,
      fixed_dt=dt, nx=nx, max_grid_size=max_grid_size)
    jobs.add(job)

  jobs.submit_all()


###############################################################################
# flameball multi-rate convergence tests

@task
def flameball_mrconv():
  """Convergence tests for the FlameBall example using MRSDC."""

  setenv('Combustion/SMC/bin/FlameBall', find_exe=True)
  rsync()
  make()

  jobs = JobQueue(rwd=env.rwd, queue='regular')

  dt0 = 5e-9
  nx  = 32

  stop_time = 40*dt0
  nprocs = 1

  #
  # reference run
  #

  dt     = dt0/4
  nnodes = 5

  name = 'nx%03d_gl%d_dt%g' % (nx, nnodes, dt)
  job  = Job(name=name, param_file='inputs-mrconv', rwd='mrconv/'+name, width=nprocs,
             walltime="01:20:00")

  job.update_params(
    nx=nx,
    fixed_dt=dt,
    stop_time=stop_time,
    cflfac=-1.0,
    sdc_nnodes=nnodes,
    sdc_nnodes_chemistry=-1,
    sdc_iters=2*nnodes-2,
    max_grid_size=nx,
    advance_method=2,
    )

  jobs.add(job)

  #
  # multi-rate runs
  #

  for nnodes, dt in product( [ (3, 9), (3, 13) ],
                             [ 2*dt0, dt0, dt0/2 ] ):

    name = 'nx%03d_gl%d.%d_dt%g' % (nx, nnodes[0], nnodes[1], dt)
    job  = Job(name=name, param_file='inputs-mrconv', rwd='mrconv/'+name, width=nprocs,
               walltime="00:40:00")

    job.update_params(
      nx=nx,
      fixed_dt=dt,
      stop_time=stop_time,
      cflfac=-1.0,
      sdc_nnodes=nnodes[0],
      sdc_nnodes_chemistry=nnodes[1],
      sdc_iters=2*nnodes[0]-2,
      max_grid_size=nx,
      advance_method=3,
      )

    #jobs.add(job)

  jobs.submit_all()
    

###############################################################################
# flamebox space/time convergence tests

@task
def flameball_stconv():
  """Convergence tests for the FlameBall example."""

  setenv('Combustion/SMC/bin/FlameBall', find_exe=True)
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


###############################################################################
# compute/plot tasks

@task
def flameball_mrconv_comp():
  """Compute convergence errors (run mrconv-comp.py) on the remote host."""
  setenv()
  rsync()

  if env.host == 'hopper' or env.host[:6] == 'edison':
    with prefix('module load numpy'):
      with cd(env.scratch + '/Combustion/SMC/analysis'):
        run('python mrconv-comp.py')
  elif env.host[:5] == 'gigan':
    with cd(env.scratch + '/Combustion/SMC/analysis'):
      run('/home/memmett/venv/base/bin/python mrconv-comp.py')


@task
def flameball_mrconv_plot():
  """Get convergence results and plot them (run mrconv-plot.py) on the local host."""
  setenv()
  get(env.scratch + 'Combustion/SMC/analysis/mrconv.pkl', 'mrconv.pkl')
  local(env.python + ' mrconv-plot.py')


@task
def flameball_stconv_comp():
  """Compute convergence errors (run stconv-comp.py) on the remote host."""
  setenv()
  rsync()

  if env.host == 'hopper' or env.host[:5] == 'edison':
    with prefix('module load numpy'):
      with cd(env.scratch + '/Combustion/SMC/analysis'):
        run('python stconv-comp.py')
  elif env.host[:5] == 'gigan':
    with cd(env.scratch + '/Combustion/SMC/analysis'):
      run('/home/memmett/venv/base/bin/python stconv-comp.py')


@task
def flameball_stconv_plot():
  """Get convergence results and plot them (run stconv-plot.py) on the local host."""
  setenv()
  get(env.scratch + 'Combustion/SMC/analysis/stconv.pkl', 'stconv.pkl')
  local(env.python + ' stconv-plot.py')


@task
def flamebox_speed_comp():
  """Compute convergence errors (run speed-comp.py) on the remote host."""
  setenv()
  rsync()

  if env.host[:6] == 'hopper' or env.host[:6] == 'edison':
    with prefix('module load numpy'):
      with cd(env.scratch + '/Combustion/SMC/analysis'):
        run('python speed-comp.py')
  elif env.host[:5] == 'gigan':
    with cd(env.scratch + '/Combustion/SMC/analysis'):
      run('/home/memmett/venv/base/bin/python speed-comp.py')


@task
def flamebox_speed_plot():
  """Get convergence results and plot them (run stconv-plot.py) on the local host."""
  setenv()
  get(env.scratch + 'Combustion/SMC/analysis/speed.pkl', 'speed.pkl')
  local(env.python + ' speed-plot.py')


###############################################################################
# build, rsync, setenv etc

@task
def flamebox_build():
  """Sync and make the FlameInABox example."""

  setenv('Combustion/SMC/bin/FlameInABox')
  make()


@task
def flameball_build():
  """Sync and make the FlameBall example."""

  setenv('Combustion/SMC/bin/FlameBall')
  make()


@task
def rsync():
  """Push (rsync) directories in env.rsync to env.host."""

  if env.host == 'localhost':
    return

  for src, dst in env.rsync:
      command = "rsync -avz -F {src}/ {host}:{dst}".format(
          host=env.host_rsync, src=src, dst=dst)
      local(command)


@task
def make(make_args=''):
  """Run make in the env.rwd directory on the remote host."""

  with cd(env.rwd):
    run('make ' + make_args)


def setenv(rwd=None, find_exe=False):
  """Setup Fabric and jobtools environment."""

  projects = '/home/memmett/projects/'

  if env.host[:6] == 'edison':
    env.scratch     = '/scratch1/scratchdirs/memmett/'
    env.scheduler   = 'edison'
    env.host_string = 'edison.nersc.gov'
    env.host_rsync  = 'edison-s'
    # env.exe         = 'main.Linux.Intel.prof.mpi.omp.exe'
    # env.exe         = 'main.Linux.Intel.debug.prof.omp.exe'
    # env.exe = 'main.Linux.gfortran.prof.mpi.omp.exe'
    env.exe = 'main.Linux.Intel.prof.mpi.omp.exe'

    env.depth   = 8
    env.pernode = 2

  elif env.host[:6] == 'hopper':
    env.scratch     = '/scratch/scratchdirs/memmett/'
    env.scheduler   = 'hopper'
    env.host_string = 'hopper.nersc.gov'
    env.host_rsync  = 'hopper-s'
    env.exe         = 'main.Linux.gfortran.prof.mpi.omp.exe'

    env.depth   = 6
    env.pernode = 4

    env.aprun_opts = [ '-cc numa_node' ]


  elif env.host[:5] == 'gigan':
    env.scratch     = '/scratch/memmett/'
    env.scheduler   = 'serial'
    env.host_string = 'gigan.lbl.gov'
    env.host_rsync  = 'gigan-s'
    env.exe         = 'main.Linux.gfortran.omp.exe'

    env.width = 1
    env.depth = 16

  if rwd:
    env.rwd = env.scratch + rwd


  if find_exe:
    with cd(env.rwd):
      exes = run('/bin/ls -1t *.exe').split()
      env.exe = exes[0]

  env.rsync = [ (projects + 'Combustion', env.scratch + 'Combustion'),
                (projects + 'BoxLib',     env.scratch + 'BoxLib'),
                (projects + 'SDCLib',     env.scratch + 'SDCLib') ]

  env.python = '/home/memmett/anaconda/bin/python'
