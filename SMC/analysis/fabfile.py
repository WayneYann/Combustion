"""Fabric (fabfile.org) tasks for SMC.  See README."""

import numpy as np

from fabric.api import *
from jobtools import JobQueue, Job
from itertools import product
from collections import namedtuple

import pickle
import compare

from pylab import *

class Container():
    pass


###############################################################################
# flamebox multi-rate speed tests

fbox = Container()

fbox.stop_time = 6e-9 * 70
fbox.dt_ref    = 1e-9
fbox.dt        = [ 6e-9, 3e-9, 1.5e-9 ]
fbox.nx        = 32
fbox.mrruns    = [ (5, 9), (7, 13) ]   # (trat, nnodes_fine)

SpeedTuple = namedtuple('SpeedTuple', [ 'scheme', 'dt', 'error', 'runtime', 'ad_evals', 'r_evals' ])


@task
def flamebox_speed():
  """Timing tests for the FlameInABox example using MRSDC."""

  setenv('Combustion/SMC/bin/FlameInABox', find_exe=True)

  jobs = JobQueue(queue='regular', walltime="01:00:00")

  max_grid_size = 32
  nprocs        = fbox.nx**3 / max_grid_size**3

  # reference run
  dt   = fbox.dt_ref
  name = 'gl3_dt%e' % dt
  job  = Job(name=name, param_file='inputs-flamebox', rwd='speed/'+name, 
            width=nprocs, walltime="02:00:00")
  job.update_params(
    advance_method=2, stop_time=fbox.stop_time,
    sdc_nnodes=3, sdc_nnodes_fine=0, sdc_iters=4,
    fixed_dt=dt, nx=fbox.nx, max_grid_size=max_grid_size).add_to(jobs)

  # timing runs
  for dt0 in fbox.dt:

    # single-rate sdc run
    dt   = dt0
    name = 'gl3_dt%e' % dt
    job = Job(name=name, param_file='inputs-flamebox', rwd='speed/'+name,
              width=nprocs, walltime="02:00:00")
    job.update_params(
      advance_method=2, stop_time=fbox.stop_time,
      sdc_nnodes=3, sdc_nnodes_fine=0, sdc_iters=4,
      fixed_dt=dt, nx=fbox.nx, max_grid_size=max_grid_size).add_to(jobs)

    # multi-rate sdc runs
    for trat, nnodes in fbox.mrruns:
      dt   = trat*dt0
      name = 'gl3.9_dt%e' % dt
      job  = Job(name=name, rwd='speed/'+name, param_file='inputs-flamebox',
                 width=nprocs)
      job.update_params(
        advance_method=3, stop_time=fbox.stop_time,
        sdc_nnodes=3, sdc_nnodes_fine=nnodes, sdc_iters=5,
        fixed_dt=dt, nx=fbox.nx, max_grid_size=max_grid_size).add_to(jobs)

  jobs.submit_all()

  
@task
def flamebox_speed_compare():
  """Compute flamebox timing and errors."""

  setenv('Combustion/SMC/bin/FlameInABox/speed')

  ref  = 'gl3_dt%e' % fbox.dt_ref
  norm = 0

  errors = []

  for dt0 in fbox.dt:

    name = 'gl3_dt%e' % dt0
    rt, nad, nr = compare.runtime(name)
    error = compare.error(name, ref, fbox.stop_time,
                          refratio=1, norm=norm, variable='pressure')

    errors.append(SpeedTuple('gl3', dt0, error, rt, nad, nr))

    for trat, nnodes in fbox.mrruns:
      dt = trat*dt0
      name = 'gl3.9_dt%e' % dt
      rt, nad, nr = compare.runtime(name)
      error = compare.error(name, ref, fbox.stop_time,
                            refratio=1, norm=norm, variable='pressure')

      errors.append(SpeedTuple('gl3.%d' % nnodes, dt0, error, rt, nad, nr))

  with open('speed.pkl', 'w') as f:
    pickle.dump(errors, f)

    
@task
def flamebox_speed_plot():
  """Plot flamebox timing and speed results."""

  with open('speed.pkl', 'r') as f:
    speed = pickle.load(f)

    schemes = set([ x.scheme for x in speed ])

    pens = {
        'gl3': { 'marker': 'o', 'ms': 12, 'color': 'black', 'label': 'SR GL 3', 'lw': 2 },
        'gl3.9': { 'marker': 's', 'ms': 12, 'color': 'blue', 'label': 'MR GL 3, 9', 'lw': 2 },
        'gl3.13': { 'marker': '^', 'ms': 12, 'color': 'red' , 'label': 'MR GL 3, 13', 'lw': 2 },
        }

    # runtime vs dt
    figure()
    for scheme in schemes:
        x, y = np.asarray(sorted([ (r.dt, r.runtime) for r in speed 
                                   if r.scheme == scheme ])).transpose()
        semilogx(x, y, **pens[scheme])

    legend()
    xlabel('dt')
    ylabel('runtime')
    savefig('mrsdc_runtime_vs_dt.eps')


    # adv/diff evals vs dt
    figure()
    for scheme in schemes:
        x, y = np.asarray(sorted([ (r.dt, r.ad_evals) for r in speed 
                                   if r.scheme == scheme ])).transpose()
        semilogx(x, y, **pens[scheme])

    legend()
    xlabel('dt')
    ylabel('number of adv/diff evals')
    savefig('mrsdc_nfevals_vs_dt.eps')

    # runtime vs error
    figure()
    for scheme in schemes:
        x, y = np.asarray(sorted([ (r.error, r.runtime) for r in speed 
                                   if r.scheme == scheme ])).transpose()
        semilogx(x, y, **pens[scheme])

    legend()
    xlabel('error')
    ylabel('runtime')
    savefig('mrsdc_error_vs_runtime.eps')


    # adv/diff evals vs error
    figure()
    for scheme in schemes:
        x, y = np.asarray(sorted([ (r.error, r.ad_evals) for r in speed 
                                   if r.scheme == scheme ])).transpose()
        semilogx(x, y, **pens[scheme])

    legend()
    xlabel('error')
    ylabel('number of adv/diff evals')
    savefig('mrsdc_nfevals_vs_error.eps')

    # error vs dt
    figure()
    for scheme in schemes:
        x, y = np.asarray(sorted([ (r.dt, r.error) for r in speed 
                                   if r.scheme == scheme ])).transpose()
        loglog(x, y, **pens[scheme])

    legend()
    xlabel('dt')
    ylabel('error')
    savefig('mrsdc_error_vs_dt.eps')
    


###############################################################################
# flameball multi-rate convergence tests

@task
def flameball_mrconv():
  """Convergence tests for the FlameBall example using MRSDC."""

  setenv('Combustion/SMC/bin/FlameBall', find_exe=True)
  rsync()

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
# build, rsync, setenv etc

@task
def build(bin):
  """Sync and build according to env.bin."""

  setenv(bin=bin)
  rsync()

  with cd(env.rwd):
    run('rm -f *.exe')
    run('make -j 4')


@task
def rsync():
  """Push (rsync) directories in env.rsync to env.host."""

  if env.host == 'localhost':
    return

  for src, dst in env.rsync:
      command = "rsync -aFvz {src}/ {host}:{dst}".format(
          host=env.host_string, src=src, dst=dst)
      local(command)


def setenv(rwd=None, bin=None, find_exe=False):
  """Setup Fabric and jobtools environment."""

  env.debug = False

  if bin is not None:
    rwd = { 'flamebox': 'Combustion/SMC/bin/FlameInABox',
            'flameball': 'Combustion/SMC/bin/FlameBall' }[bin]

  projects = '/home/memmett/projects/'

  if env.host[:6] == 'edison':
    env.scratch     = '/scratch1/scratchdirs/memmett/'
    env.scheduler   = 'edison'
    env.host_string = 'edison.nersc.gov'

    env.depth   = 8
    env.pernode = 2

  elif env.host[:6] == 'hopper':
    env.scratch     = '/scratch/scratchdirs/memmett/'
    env.scheduler   = 'hopper'
    env.host_string = 'hopper.nersc.gov'

    env.depth   = 6
    env.pernode = 4

    env.aprun_opts = [ '-cc numa_node' ]

  elif env.host[:5] == 'gigan':
    env.scratch     = '/scratch1/memmett/'
    env.scheduler   = 'serial'
    env.host_string = 'gigan.lbl.gov'
    env.ffdcompare  = env.scratch + 'AmrPostprocessing/F_Src/ffdcompare.Linux.gfortran.exe'

    env.width = 1
    env.depth = 8

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
