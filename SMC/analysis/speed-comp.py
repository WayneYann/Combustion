"""Compute space-time errors."""

from __future__ import with_statement

import pickle
import compare

compare.setenv()


def speed_comp():

  compare.env.rwd = 'Combustion/SMC/bin/FlameInABox/speed'

  dt0 = 6e-9
  ref = 'gl3_dt1.000000e-09'
  time = 4e-8 + dt0*70

  errors = {}

  norm = 0

  for dt0 in [ 6e-9, 3e-9, 1.5e-9 ]:

    dt  = dt0
    run = 'gl3_dt%e' % dt

    timing = compare.runtime(run)
    error  = compare.error(run, ref, time, refratio=1, norm=norm, variable='pressure')

    errors['gl3', dt] = error, timing

    dt  = 5*dt0
    run = 'gl3.9_dt%e' % dt

    timing = compare.runtime(run)
    error  = compare.error(run, ref, time, refratio=1, norm=norm, variable='pressure')

    errors['gl3.9', dt] = error, timing

    dt  = 7*dt0
    run = 'gl3.13_dt%e' % dt

    timing = compare.runtime(run)
    error  = compare.error(run, ref, time, refratio=1, norm=norm, variable='pressure')

    errors['gl3.13', dt] = error, timing


  with open('speed.pkl', 'w') as f:
    pickle.dump(errors, f)


if __name__ == '__main__':
  speed_comp()
