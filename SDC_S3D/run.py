"""SDC_S3D run script."""

from runutils import autorun

defaults = {
  'tfinal':  2.5e-7, 
  'dt':      2.5e-7/10, 
  'cfl':    -1.0,             # fixed step size
  'method':  2,
  'nnodes':  3,
  }

autorun('test', **defaults)


