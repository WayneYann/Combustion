
from fabric.api import *

import os
import re
from StringIO import StringIO

def memoize(obj):
    import functools
    
    cache = obj.cache = {}

    @functools.wraps(obj)
    def memoizer(*args, **kwargs):
        key = str(args) + str(kwargs)
        if key not in cache:
            cache[key] = obj(*args, **kwargs)
        return cache[key]
    return memoizer

  
@memoize
def find_plotfile(rundir, time):
  """Return plotfile with time *time* in run directory *rundir*."""

  with cd(os.path.join(env.rwd, rundir)):
    program = [ "function abs(x) { if (x<0.0) return -x; return x }",
                "FNR==2 { ncomp=$1 }",
                "FNR==ncomp+4 { if (abs($1-%e) < 1e-12) print(FILENAME) }" % time ]
    header = run("awk '%s' plt*/Header" % '; '.join(program), quiet=True)
    if header:
      return os.path.join(rundir, header.split('/')[0])

  return None


def compare(p1, p2, variables, norm=0):

  with cd(env.rwd):
    output = run('%s --infile1 %s --infile2 %s' % (env.ffdcompare, p1, p2), quiet=True)
  
  err = {}
  for l in output.split('\n'):
    r = l.split()
    if r and r[0] in variables:
      err[r[0]] = float(r[norm+1])

  return err


def runtime(rundir):

  with cd(os.path.join(env.rwd, rundir)):
    out = run('grep -e "Total Run Time =" -e "AD fevals =" -e "R  fevals =" stdout', quiet=True)

  m = re.search(r"Total Run Time =\s*(\S*)", out)
  time = float(m.group(1)) if m else -1.0

  m = re.search(r"AD fevals =\s*(\S*)", out)
  nad = int(m.group(1)) if m else -1

  m = re.search(r"R  fevals =\s*(\S*)", out)
  nr = int(m.group(1)) if m else -1

  return time, nad, nr


def error(rundir1, rundir2, time, norm=2, variables=['pressure'], refratio=1, diff=None):
  """Compute the error of *variable* between runs at time *time*."""

  p1 = find_plotfile(rundir1, time)
  p2 = find_plotfile(rundir2, time)

  if p1 is None or p2 is None:
    print '  plotfiles not found!'
    return None

  err = compare(p1, p2, variables, norm)

  if len(variables) == 1:
    return err[variables[0]]

  return err
