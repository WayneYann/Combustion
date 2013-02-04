"""Fabric (fabfile.org) utilities for submitting jobs."""

import os

from fabric.api import *
from fabric.colors import *
from fabric.utils import *
from fabric.contrib.files import *

from textwrap import dedent


def submit(*args, **kwargs):

  import sys
  thismodule = sys.modules[__name__]

  sub = getattr(thismodule, 'submit_' + env.host)
  sub(*args, **kwargs)


def submit_hopper(name=None, rundir=None, probin=None, queue=None,
                  nprocs=None, nthreads=None, pernode=None,
                  walltime='00:30:00', **kwargs):

  if nprocs is None:
    nprocs = env.nprocs

  if nthreads is None:
    nthreads = env.nthreads

  if pernode is None:
    pernode = env.pernode

  if queue is None:
    queue = 'regular'

  pbs = '{bin}/{rundir}/{name}.pbs'.format(bin=env.bin, rundir=rundir, name=name)

  puts(green('generating: ' + name))
  run(dedent("""\
             cat > {pbs} << EOF
             #!/bin/bash -l
             #PBS -N {name}
             #PBS -q {queue}
             #PBS -l mppwidth={nprocs}
             #PBS -l mppnppn={pernode}
             #PBS -l mppdepth={nthreads}
             #PBS -l walltime={walltime}
             #PBS -o {bin}/{rundir}/out
             #PBS -e {bin}/{rundir}/err
             #PBS -V

             cd {bin}/{rundir}

             export OMP_NUM_THREADS={nthreads}

             aprun -B {bin}/{exe} {probin}
             EOF
             """.format(
               name      = name,
               host      = env.host,
               bin       = env.bin,
               exe       = env.exe,
               rundir    = rundir,
               nprocs    = nprocs,
               nthreads  = nthreads,
               pernode   = pernode,
               queue     = queue,
               walltime  = walltime,
               probin    = probin,
               pbs       = pbs)))

  puts(green('submitting: ' + name))
  run('qsub ' + pbs)


def submit_gigan(name=None, rundir=None, probin=None, nprocs=None, 
                 nthreads=None, **kwargs):

  if nprocs is None:
    nprocs = env.nprocs

  if nthreads is None:
    nthreads = env.nthreads

  pbs = '{bin}/{rundir}/{name}.sh'.format(bin=env.bin, rundir=rundir, name=name)

  puts(green('generating: ' + name))
  run(dedent("""\
             cat > {pbs} << EOF
             #!/bin/bash -l
             #PBS -N {name}
             #PBS -q {queue}
             #PBS -l mppwidth={nprocs}
             #PBS -l mppnppn={pernode}
             #PBS -l mppdepth={nthreads}
             #PBS -l walltime={walltime}
             #PBS -o {bin}/{rundir}/out
             #PBS -e {bin}/{rundir}/err
             #PBS -V

             cd {bin}/{rundir}

             export OMP_NUM_THREADS={nthreads}

             mpirun -n {nprocs} {bin}/{exe} {probin}
             EOF
             """.format(
               name      = name,
               host      = env.host,
               bin       = env.bin,
               exe       = env.exe,
               rundir    = rundir,
               nprocs    = nprocs,
               nthreads  = nthreads,
               probin    = probin,
               pbs       = pbs)))

  puts(green('submitting: ' + name))
  run('sh ' + pbs)
