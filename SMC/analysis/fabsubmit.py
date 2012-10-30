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


def submit_hopper(name, rundir, probin, queue=None, walltime='00:30:00', **kwargs):

  if not queue:
    queue = 'regular'

  puts(green('qsub script: ' + name))
  run(dedent("""\
             cat > {name}.pbs << EOF
             #!/bin/bash -l
             #PBS -N {rundir}
             #PBS -q {queue}
             #PBS -l mppwidth={nproc}
             #PBS -l mppnppn={pernode}
             #PBS -l mppdepth={nthreads}
             #PBS -l walltime={walltime}
             #PBS -k oe
             #PBS -V

             cd {bin}/{rundir}

             export OMP_NUM_THREADS={nthreads}

             aprun -B {bin}/{exe} {probin}
             EOF
             """.format(
               host      = env.host,
               bin       = env.bin,
               exe       = env.exe,
               rundir    = rundir,
               nproc     = env.nprocs,
               nthreads  = env.nthreads,
               pernode   = env.pernode,
               name      = name,
               queue     = queue,
               walltime  = walltime,
               probin    = probin)))

  puts(green('submitting: ' + name))
  run('qsub {name}.pbs'.format(name=name))
