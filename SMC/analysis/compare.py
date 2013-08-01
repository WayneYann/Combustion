
from __future__ import with_statement

import glob
import os
import re
import subprocess


class Container(object):
    pass
env = Container()

def setenv(rwd=None, debug=False):
    import socket
    host = socket.gethostname()

    env.debug = debug

    if host == 'gigan':
        env.scratch    = '/scratch/memmett/Combustion/SMC/bin/FlameBall/'
        env.ffdcompare = '/scratch/memmett/AmrPostprocessing/F_Src/ffdcompare.Linux.gfortran.exe'
    elif host[:6] == 'hopper':
        env.scratch    = '/scratch/scratchdirs/memmett/'
        env.ffdcompare = '/global/homes/m/memmett/projects/AmrPostprocessing/F_Src/ffdcompare.Linux.gfortran.exe'

    elif host[:6] == 'edison':
        env.scratch    = '/scratch1/scratchdirs/memmett/'
        env.ffdcompare = '/global/homes/m/memmett/projects/AmrPostprocessing/F_Src/ffdcompare.Linux.gfortran.exe'

    if rwd:
        env.rwd = rwd



def find_plotfile(rundir, time):
  """Return plotfile with time *time* in run directory *rundir*."""

  plts = glob.glob(os.path.join(env.scratch, env.rwd, rundir, 'plt*'))
  if env.debug:
      print rundir, len(plts)
  for plt in plts:
    with open(plt + '/Header', 'r') as f:
      header = f.read().split('\n')
    ncomp = int(header[1])
    plttime = float(header[ncomp+3])
    if env.debug:
        print rundir, plttime, time
    if abs(plttime - time) < 1e-18:
      return plt

  return None


def compare(p1, p2, variable, norm=0):

  proc = subprocess.Popen([env.ffdcompare, '--infile1', p1, '--infile2', p2], stdout=subprocess.PIPE)
  output, _ = proc.communicate()
  
  for l in output.split('\n'):
    r = l.split()
    if r and r[0] == variable:
      return float(r[norm+1])

  return None


def runtime(run):

  with open(os.path.join(env.scratch, env.rwd, run, 'stdout'), 'r') as f:
    stdout = f.read()

  m = re.search(r"Total Run Time =\s*(\S*)", stdout)
  time = float(m.group(1)) if m else -1.0

  m = re.search(r"AD fevals =\s*(\S*)", stdout)
  nad = int(m.group(1)) if m else -1

  m = re.search(r"R  fevals =\s*(\S*)", stdout)
  nr = int(m.group(1)) if m else -1

  return time, nad, nr


def error(rundir1, rundir2, time, norm=2, variable='pressure', refratio=1, diff=None):
  """Compute the error of *variable* between runs at time *time*."""


  print 'computing errors:', rundir1, rundir2, time, variable, refratio

  p1 = find_plotfile(rundir1, time)
  p2 = find_plotfile(rundir2, time)

  if p1 is None or p2 is None:
    print '  plotfiles not found!'
    return None

  err = compare(p1, p2, variable, norm)

  print '  p1:', p1
  print '  p2:', p2
  print '  l%d:' % norm, err

  return err
