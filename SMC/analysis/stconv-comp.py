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
    if abs(plttime - time) < 1e-8:
      return plt

  return None


def error(rundir1, rundir2, time, variable='pressure', refratio=1, diff=None):
  """Compute the error of *variable* between runs at time *time*."""

  print 'computing errors:', rundir1, rundir2, time, variable, refratio

  p1 = find_plotfile(rundir1, time)
  p2 = find_plotfile(rundir2, time)

  if p1 is None or p2 is None:
    print '  plotfiles not found!'
    return None

  errs, dnames = compare(p1, p2, refratio=refratio, variables=[variable], diff=diff, norm=0)

  print '  p1:', p1
  print '  p2:', p2
  print '  l2:', errs[variable][0]

  return errs[variable][0]


def flameball_stconv_comp():

  import stconv

  env.base = '/scratch/scratchdirs/memmett/Combustion/SMC/bin/FlameBall/stconv2'

  errors = { 'stconv': {}, 'cflconv': {} }
  for nx, dt, nnodes in stconv.runs:

    # compare to same grid run
    refdir   = 'nx%03d_gl%d_dt%g' % (nx, 5, stconv.min_dt[nx])
    rundir   = 'nx%03d_gl%d_dt%g' % (nx, nnodes, dt)

    if refdir != rundir:
      err = error(rundir, refdir, stconv.stop_time, refratio=1, variable='pressure')
      if err:
        errors['cflconv'][nx, dt, nnodes] = err

    # compare to fine grid run
    refdir   = 'nx%03d_gl%d_dt%g' % (128, 5, stconv.min_dt[128])
    rundir   = 'nx%03d_gl%d_dt%g' % (nx, nnodes, dt)

    if refdir != rundir:
      err = error(rundir, refdir, stconv.stop_time, refratio=128/nx, variable='pressure', 
                  diff='diff_nx%03d_gl%d_dt%g_box' % (nx, nnodes, dt))
      if err:
        errors['stconv'][nx, dt, nnodes] = err

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
