"""Compute space-time errors."""

from __future__ import with_statement
from itertools import product

import pickle
import compare

compare.setenv('Combustion/SMC/bin/FlameBall', False)


def flameball_mrconv_comp():

  errors = { 'dtconv': {} }

  dt0  = 5e-9
  nx   = 32
  time = 40 * dt0

  for nodes, dt in product( [ (3, 9), (3, 13) ],
                             [ 2*dt0, dt0, dt0/2 ] ):

      ref = 'mrconv/nx%03d_%s_dt%g' % (nx, 'gl5', dt0/4)
      run = 'mrconv/nx%03d_gl%d.%d_dt%g' % (nx, nodes[0], nodes[1], dt)

      try:
          timing = compare.runtime(run)
          error  = compare.error(run, ref, time, refratio=1, norm=2, variable='pressure')
      except IOError:
          timing = (0.0, 0, 0)
          error  = 0.0
          
      if error:
          errors['dtconv'][nx, dt, nodes] = error, timing

  with open('mrconv.pkl', 'w') as f:
    pickle.dump(errors, f)


if __name__ == '__main__':
  flameball_mrconv_comp()
