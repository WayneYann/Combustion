"""Compute space-time errors."""

import glob
import os
import re
import pickle

import numpy as np


from itertools import product
from collections import namedtuple
from fdcompare import fdcompare as compare

class Container():
  pass

env = Container()


def find_plotfile(rundir, time):
  """Return plotfile with time *time* in run directory *rundir*."""

  plts = glob.glob(os.path.join(env.base, rundir, 'plt*'))
  for plt in plts:
    with open(plt + '/Header', 'r') as f:
      header = f.read().split('\n')
    ncomp = int(header[1])
    plttime = float(header[ncomp+3])
    if abs(plttime - time) < 1e-13:
      return plt

  return None


def error(rundir1, rundir2, time, variable='pressure', refratio=1):
  """Compute the error of *variable* between runs at time *time*."""

  print 'computing errors:', rundir1, rundir2, time, variable, refratio

  p1 = find_plotfile(rundir1, time)
  p2 = find_plotfile(rundir2, time)

  if p1 is None or p2 is None:
    print '  plotfiles not found!'
    return None

  errs, dnames = compare(p1, p2, refratio=refratio, variables=[variable])

  print '  p1:', p1
  print '  p2:', p2
  print '  l2:', errs[variable][0]

  return errs[variable][0]


def flameball_stconv_comp():

  env.base = '/scratch/scratchdirs/memmett/Combustion/SMC/bin/FlameBall/stconv'

  dt0 = 1e-7
  stop_time = dt0 * 10

  NX     = [ 32, 64, 128 ]
  DT     = [ dt0/8, dt0/4, dt0/2, dt0 ]
  NNODES = [ 3, 5 ]

  errors = { 'stconv': {}, 'cflconv': {} }
  for nx, dt, nnodes in product(NX, DT, NNODES):

    # compare to same grid run
    refdir   = 'nx%03d_gl%d_dt%g' % (nx, max(NNODES), min(DT))
    rundir   = 'nx%03d_gl%d_dt%g' % (nx, nnodes, dt)

    if refdir != rundir:
      err = error(rundir, refdir, stop_time, refratio=1, variable='density')
      if err:
        errors['cflconv'][nx, dt, nnodes] = err

  print errors
  with open('stconv.pkl', 'w') as f:
    pickle.dump(errors, f)



def flameball_cflconv_comp(nx=32):

  env.base = '/scratch/scratchdirs/memmett/Combustion/SMC/bin/FlameBall/stconv'

  errors = {}
  for cfl, nnodes in product( [ 0.5, 0.75, 1.0 ], [ 3, 5 ] ):

    refdir = 'nx%03d_gl5_cfl0.25' % nx
    rundir = 'nx%03d_gl%d_cfl%.2f' % (nx, nnodes, cfl)

    print 'computing errors:', rundir
    errors[nx, cfl, nnodes] = error(rundir, refdir)


  with open('cflconv%03d.pkl' % nx, 'w') as f:
    pickle.dump(errors, f)



if __name__ == '__main__':
  flameball_stconv_comp()
  #flameball_cflconv_comp(32)
  #flameball_cflconv_comp(64)
