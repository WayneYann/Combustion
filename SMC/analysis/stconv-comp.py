"""Compute space-time errors."""

import pickle
import glob
import os
import re

from itertools import product
from collections import namedtuple

import numpy as np
import sh


class Container():
  pass

env = Container()
key = namedtuple('key', [ 'nx', 'cfl', 'nodes' ])

comp_dname = '/global/homes/m/memmett/projects/BoxLib/Tools/C_util/Convergence/'
comp_fname = 'DiffSameDomainRefined3d.Linux.Cray.Cray.ex'
comp = sh.Command(comp_dname + comp_fname)



def last_plotfile(rundir):
  """Return last plotfile in run directory."""
  return sorted(glob.glob(os.path.join(env.base, rundir, 'plt*')))[-1]


def error(rundir1, rundir2):

  p1  = last_plotfile(rundir1)
  p2  = last_plotfile(rundir2)
  out = comp("infile1=%s" % p1, "reffile=%s" % p2)

  # parse output and grab errors
  m = re.search(r'^  0 (.*)$', str(out), re.MULTILINE)
  try:
    errs = m.group(1).split()
  except:
    errs = [ np.nan ]

  return max(errs)


def flameball_stconv_comp():

  env.base = '/scratch/scratchdirs/memmett/Combustion/SMC/bin/FlameBall/stconv'

  errors = {}
  for nx, cfl, nnodes in product( [ 32, 64, 128 ],
                                  [ 0.25, 0.5, 0.75, 1.0 ],
                                  [ 3, 5 ] ):

    refdir = 'nx128_gl5_cfl0.25'
    rundir = 'nx%03d_gl%d_cfl%.2f' % (nx, nnodes, cfl)

    print 'computing errors:', rundir
    errors[key(nx, cfl, nnodes)] = error(rundir, refdir)

  with open('stconv.pkl', 'w') as f:
    pickle.dump(errors, f)



if __name__ == '__main__':
  flameball_stconv_comp()
