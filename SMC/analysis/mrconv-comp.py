"""Compute space-time errors."""

from __future__ import with_statement

import pickle
import compare

compare.setenv()


def flameball_mrconv_comp():

  dt0  = 5e-9

  errors = { 'dtconv': {} }

  nx     = 32
  time   = dt0 * 10

  for nodes in [ 'rk3', 'gl3', 'gl3-9' ]:
    for dt in [ dt0, dt0/2, dt0/4, dt0/8 ]:
      ref = 'mrconv/nx%03d_%s_dt%g' % (nx, 'gl3', dt0/16)
      run = 'mrconv/nx%03d_%s_dt%g' % (nx, nodes, dt)

      try:
          timing = compare.runtime(run)
          error  = compare.error(run, ref, time, refratio=1, norm=2, variable='density')
      except IOError:
          timing = (0.0, 0, 0)
          error  = 0.0
          
      if error:
          errors['dtconv'][nx, dt, nodes] = error, timing

  with open('mrconv.pkl', 'w') as f:
    pickle.dump(errors, f)


if __name__ == '__main__':
  flameball_mrconv_comp()
