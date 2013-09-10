"""Fabric (fabfile.org) tasks for SMC.  See README."""

import numpy as np

from fabric.api import *
from fabric.colors import *
from jobtools import JobQueue, Job
from itertools import product
from collections import namedtuple
from pprint import pprint

import pickle
import compare

from pylab import *

class Container():
    pass


###############################################################################
# flamebox multi-rate speed tests

fbox = Container()

fbox.stop_time = 6e-9 * 70
fbox.dt_ref    = 1.5e-9
fbox.dt        = [ 6e-9, 3e-9, 1.5e-9 ]
fbox.nx        = 32
fbox.mrruns    = [ (5, 9), (7, 13) ]
fbox.mrreps    = [ (5, 3, 4), (7, 3, 6), (10, 3, 8), (10, 5, 2) ]

SpeedTuple = namedtuple('SpeedTuple',
                        [ 'scheme', 'dt', 'error', 'runtime', 'ad_evals', 'r_evals' ])

pens = {
    'gl3':     { 'marker': 'o', 'ms': 12, 'color': 'black', 'label': 'SR GL 3', 'lw': 2 },
    'gl3.3r4': { 'marker': '^', 'ms': 12, 'color': 'black', 'label': 'SR GL 3x4', 'lw': 2 },
    'gl3.3r6': { 'marker': 'v', 'ms': 12, 'color': 'black', 'label': 'SR GL 3x6', 'lw': 2 },
    'gl3.3r8': { 'marker': 's', 'ms': 12, 'color': 'black', 'label': 'SR GL 3x8', 'lw': 2 },
    'gl3.5r2': { 'marker': 'd', 'ms': 12, 'color': 'black', 'label': 'SR GL 5x2', 'lw': 2 },
    'gl3.9':   { 'marker': 's', 'ms': 12, 'color': 'blue',  'label': 'MR GL 3, 9', 'lw': 2 },
    'gl3.13':  { 'marker': '^', 'ms': 12, 'color': 'red' ,  'label': 'MR GL 3, 13', 'lw': 2 },
    }

@task
def flamebox_speed():
  """Timing tests for the FlameInABox example using MRSDC."""

  setenv('Combustion/SMC/bin/FlameInABox', find_exe=True)

  jobs = JobQueue(queue='regular', walltime="01:00:00")

  max_grid_size = 16
  nprocs        = fbox.nx**3 / max_grid_size**3

  # reference run
  dt   = fbox.dt_ref
  name = 'gl3_dt%e' % dt
  job  = Job(name=name, param_file='inputs-flamebox', rwd='speed/'+name, 
            width=nprocs, walltime="02:00:00")
  job.update_params(
    advance_method=2, stop_time=fbox.stop_time,
    sdc_nnodes=3, sdc_nnodes_fine=0, sdc_iters=4, repeat=1,
    fixed_dt=dt, nx=fbox.nx, max_grid_size=max_grid_size)#.add_to(jobs)

  # timing runs
  for dt0 in fbox.dt:

    # single-rate sdc run
    dt   = dt0
    name = 'gl3_dt%e' % dt
    job = Job(name=name, param_file='inputs-flamebox', rwd='speed/'+name,
              width=nprocs, walltime="02:00:00")
    job.update_params(
      advance_method=2, stop_time=fbox.stop_time,
      sdc_nnodes=3, sdc_nnodes_fine=0, sdc_iters=4, repeat=1,
      fixed_dt=dt, nx=fbox.nx, max_grid_size=max_grid_size).add_to(jobs)

    # multi-rate sdc runs
    for trat, nnodes in fbox.mrruns:
      dt   = trat*dt0
      name = 'gl3.%d_dt%e' % (nnodes, dt)
      job  = Job(name=name, rwd='speed/'+name, param_file='inputs-flamebox',
                 width=nprocs)
      job.update_params(
        advance_method=3, stop_time=fbox.stop_time,
        sdc_nnodes=3, sdc_nnodes_fine=nnodes, sdc_iters=4, repeat=1,
        fixed_dt=dt, nx=fbox.nx, max_grid_size=max_grid_size).add_to(jobs)

    for trat, nnodes, nrep in fbox.mrreps:
      dt   = trat*dt0
      name = 'gl3.%dr%d_dt%e' % (nnodes, nrep, dt)
      job  = Job(name=name, rwd='speed/'+name, param_file='inputs-flamebox',
                 width=nprocs)
      job.update_params(
        advance_method=3, stop_time=fbox.stop_time,
        sdc_nnodes=3, sdc_nnodes_fine=nnodes, sdc_iters=4, repeat=nrep,
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
                          refratio=1, norm=norm, variables=['pressure'])
    errors.append(SpeedTuple('gl3', dt0, error, rt, nad, nr))

    for trat, nnodes in fbox.mrruns:
      dt = trat*dt0
      name = 'gl3.%d_dt%e' % (nnodes, dt)
      rt, nad, nr = compare.runtime(name)
      error = compare.error(name, ref, fbox.stop_time,
                            refratio=1, norm=norm, variables=['pressure'])
      errors.append(SpeedTuple('gl3.%d' % nnodes, dt, error, rt, nad, nr))

    for trat, nnodes, nrep in fbox.mrreps:
      dt   = trat*dt0
      name = 'gl3.%dr%d_dt%e' % (nnodes, nrep, dt)
      rt, nad, nr = compare.runtime(name)
      error = compare.error(name, ref, fbox.stop_time,
                            refratio=1, norm=norm, variables=['pressure'])
      errors.append(SpeedTuple('gl3.%dr%d' % (nnodes, nrep), dt, error, rt, nad, nr))

  with open('speed.pkl', 'w') as f:
    pickle.dump(errors, f)

    
@task
def flamebox_speed_plot(interactive=False):
  """Plot flamebox timing and speed results."""

  with open('speed.pkl', 'r') as f:
    speed = pickle.load(f)

  schemes = set([ x.scheme for x in speed ])

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

  pprint(sorted([ (r.scheme, r.dt, r.r_evals, r.ad_evals, r.runtime, r.error) for r in speed ]))

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

  if interactive:
    show()

###############################################################################
# flameball multi-rate convergence tests

fball = Container()

dt = 5e-9

fball.stop_time = dt * 40
fball.dt_ref    = dt / 4
fball.dt        = [ 2*dt, dt, dt/2 ]
fball.nx        = 32
fball.mrruns    = [ (3, 9, 1), (3, 13, 1), (3, 5, 2) ]
fball.variables = [
    'density',
    'temperature',
    'x_vel',
    # 'y_vel',
    # 'z_vel',
    'Y(H2)',
    'Y(O2)',
    'Y(OH)',
    'Y(H2O)',
    'Y(N2)',
    ]

texvar = {
    'density': r'$\rho$',
    'temperature': r'$T$',
    'x_vel': r'$u_x$',
    'y_vel': r'$u_y$',
    'z_vel': r'$u_z$',
    'Y(H2)': r'$Y({\rm H}_2)$',
    'Y(O2)': r'$Y({\rm O}_2)$',
    'Y(OH)': r'$Y({\rm OH})$',
    'Y(H2O)': r'$Y({\rm H}_2{\rm O})$',
    'Y(N2)': r'$Y({\rm N}_2)$',
}


@task
def flameball_mrconv():
  """Convergence tests for the FlameBall example using MRSDC."""

  setenv('Combustion/SMC/bin/FlameBall', find_exe=True)

  max_grid_size = 16
  nprocs        = fbox.nx**3 / max_grid_size**3

  jobs = JobQueue(queue='regular', walltime="00:10:00", width=nprocs)

  # reference run
  dt   = fball.dt_ref
  name = 'gl5_dt%e' % dt
  job  = Job(name=name, param_file='inputs-mrconv', rwd='mrconv/'+name)
  job.update_params(
    advance_method=2, stop_time=fball.stop_time,
    sdc_nnodes=5, sdc_nnodes_fine=-1, sdc_iters=8, repeat=1,
    fixed_dt=dt, nx=fball.nx, max_grid_size=max_grid_size)
  jobs.add(job)

  # multi-rate runs
  for dt in fball.dt:
    for coarse, fine, repeat in fball.mrruns:
      name = 'gl%d.%dr%d_dt%e' % (coarse, fine, repeat, dt)
      job  = Job(name=name, param_file='inputs-mrconv', rwd='mrconv/'+name)
      job.update_params(
        advance_method=2, stop_time=fball.stop_time,
        sdc_nnodes=coarse, sdc_nnodes_fine=fine, sdc_iters=2*coarse-2, repeat=repeat,
        fixed_dt=dt, nx=fball.nx, max_grid_size=max_grid_size)
      jobs.add(job)

  jobs.submit_all()

@task
def flameball_mrconv_compare():
  """Compute flameball errors."""

  setenv('Combustion/SMC/bin/FlameBall/mrconv')

  ref  = 'gl5_dt%e' % fball.dt_ref

  errors = {}
  for dt in fball.dt:
    for coarse, fine, repeat in fball.mrruns:
      name = 'gl%d.%dr%d_dt%e' % (coarse, fine, repeat, dt)
      errors[dt,name,0] = compare.error(name, ref, fball.stop_time,
                                        refratio=1, norm=0, variables=fball.variables)
      errors[dt,name,2] = compare.error(name, ref, fball.stop_time,
                                        refratio=1, norm=2, variables=fball.variables)
  with open('mrconv.pkl', 'w') as f:
    pickle.dump(errors, f)

@task
def flameball_mrconv_tabulate():
    
  with open('mrconv.pkl', 'r') as f:
    errors = pickle.load(f)

  for coarse, fine, repeat in fball.mrruns:
    for i, variable in enumerate(fball.variables):
      for k, dt in enumerate(fball.dt):
        name = 'gl%d.%dr%d_dt%e' % (coarse, fine, repeat, dt)
        l0 = errors[dt,name,0][variable]
        l2 = errors[dt,name,2][variable]

        if i == 0 and k == 0:
          row = 'GL %d / GL %d $\\times$n %d & ' % (coarse, fine, repeat)
        else:
          row = ' & '
        
        if k == 0:
          row += '%s & %.2f & %.3e & & %.3e & \\\\' % (texvar[variable], dt/1e-9, l0, l2)
        else:
          name0 = 'gl%d.%dr%d_dt%e' % (coarse, fine, repeat, fball.dt[k-1])
          e   = fball.dt[k-1] / dt
          l0r = np.log(errors[fball.dt[k-1],name0,0][variable] / l0) / np.log(e)
          l2r = np.log(errors[fball.dt[k-1],name0,2][variable] / l2) / np.log(e)

          row += ' & %.2f & %.3e & %.2f & %.3e & %.2f \\\\' % (dt/1e-9, l0, l0r, l2, l2r)

        print row


###############################################################################
# build, rsync, setenv etc

@task
def build(bin, opts=''):
  """Sync and build according to env.bin."""

  setenv(bin=bin)
  rsync()

  with cd(env.rwd):
    run('rm -f *.exe')
    run('make -j 4 %s' % opts)


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
    env.ffdcompare  = env.scratch + 'AmrPostprocessing/F_Src/ffdcompare.Linux.gfortran.exe'

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


@task
def echo_env():
  setenv()
  print green("=== remote environment variables ===")
  run('env')
  print green("=== remote modules ===")
  run('module list')
  
