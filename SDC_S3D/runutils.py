"""SDC_S3D run utilities."""

import os
import subprocess

dirname    = os.path.dirname(os.path.realpath(__file__)) + '/'
executable = dirname + 'main.Linux.gfortran.mpi.omp.exe'


###############################################################################

def autorun(name, defaults={}, **arargs):

  clargs = parse_arguments()

  # make directory and change into it
  try:
    os.system('mkdir -p ' + name)
  except:
    pass

  # write probin
  probin = {}
  probin.update(defaults)
  probin.update(arargs)
  probin.update(clargs)
  
  namelist(name + '/input.nml', 'probin', **probin)

  # run!
  p = subprocess.Popen([ executable, 'input.nml' ], cwd=dirname+name)
  p.wait()
  

###############################################################################

def namelist(filename, namelist, **kwargs):
  """Write a Fortran namelist."""

  with open(filename, 'w') as f:
    f.write('&%s\n' % namelist)

    for kw in kwargs:
      f.write('  %s = %s\n' % (kw, kwargs[kw]))
    
    f.write('/\n')


###############################################################################

def parse_arguments(argv=None):
  """Return dictionary of keyword arguments extracted from arguments."""

  if argv is None:
    import sys
    argv = sys.argv

  kwargs = {}
  for a in argv:
    if a.count('=') < 1:
      continue

    k, v = a.split('=', 1)

    # integer?
    try:
      kwargs[k] = int(v)
      continue
    except:
      pass

    # float?
    try:
      kwargs[k] = float(v)
      continue
    except:
      pass

    # boolean?
    if v.lower() == 'true':
      kwargs[k] = True
      continue

    if v.lower() == 'false':
      kwargs[k] = False
      continue

    # eval?
    try:
      kwargs[k] = eval(v)
      continue
    except:
      pass

    # default to string
    kwargs[k] = v

  return kwargs




