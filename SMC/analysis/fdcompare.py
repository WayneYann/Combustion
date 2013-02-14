"""3D finite difference plotfile comparison routines."""

from pyboxlib.plotfile import plotfile
import numpy as np


###############################################################################
# helpers

def coarse_index(idx, refratio, pbxrange):
  """Compute the coarse index corresponding to the fine index *idx*
  given the refinement ration *refratio*.

  A status is also returned, and will be False if the fine index is
  not divisible by the refinement ratio or the resulting coarse index
  lies outside of the box *pbxrange*.
  """

  if idx % refratio:
    return -1, False
  if idx / refratio in pbxrange:
    return idx / refratio - pbxrange[0], True
  return -1, False


def fdcompare_boxes(plt1, comp1, plt2, comp2, level, box, norm=2):
  """Compare components between two plotfiles that have the same
  spatial resolution."""

  plt1.bind(level, box, comp1)
  plt2.bind(level, box, comp2)

  fab1 = plt1.fab(level, box)
  fab2 = plt2.fab(level, box)

  diff  = fab1.array - fab2.array
  if norm:
    error = np.sum(abs(diff)**norm)
  else:
    error = np.amax(abs(diff))

  plt2.unbind(level, box)
  plt1.unbind(level, box)
          
  return error


def fdcompare_boxes_refined(plt1, comp1, plt2, comp2, level, box1, refratio, norm, save_diff):
  """Compare components between two plotfiles that have different
  spatial resolutions.

  Note that plt1, comp1, and box1 correspond to the coarse plotfile.
  """

  plt1.bind(level, box1, comp1)
  fab1 = plt1.fab(level, box1)

  fab2coarse = -6666.6 * np.ones(fab1.shape)

  # cycle through all fine boxes and look for indexes that match the
  # coarse box we're in
  for box2 in range(1, plt2.nboxes(level)+1):

    plt2.bind(level, box2, comp2)
    fab2 = plt2.fab(level, box2)

    for i in fab2.pbxrange(1):
      i2, valid = coarse_index(i, refratio, fab1.pbxrange(1))
      if not valid: continue

      for j in fab2.pbxrange(2):
        j2, valid = coarse_index(j, refratio, fab1.pbxrange(2))
        if not valid: continue

        for k in fab2.pbxrange(3):
          k2, valid = coarse_index(k, refratio, fab1.pbxrange(3))
          if not valid: continue

          fab2coarse[i2,j2,k2] = fab2[i,j,k]

    plt2.unbind(level, box2)

  if (np.any(fab2coarse == -6666.6)):
    raise ValueError('BAD FILL.')

  diff  = fab1.array - fab2coarse
  if norm:
    error = np.sum(abs(diff)**norm) / diff.size
  else:
    error = np.amax(abs(diff))

  plt1.unbind(level, box1)

  if save_diff:
    fname = save_diff + str(box1) + '.npy'
    np.save(fname, diff)

  return error


def fdcompare(dname1, dname2, norm=2, refratio=1, variables=None, diff=None):
  """Compare two 3D finite differenced plotfiles."""

  plt1 = plotfile()
  plt1.create(dname1)

  plt2 = plotfile()
  plt2.create(dname2)

  if plt1.dim != plt2.dim:
    raise ValueError("Number of dimensions don't match: %s, %s" % (dname1, dname2))

  if plt1.flevel != plt2.flevel:
    raise ValueError("Number of levels don't match: %s, %s" % (dname1, dname2))

  if not variables:
    variables = sorted(plt1.variables)

  errors = {}
  
  for variable in variables:

    errors[variable] = {}

    comp1 = plt1.variables[variable]
    try:
      comp2 = plt2.variables[variable]
    except KeyError:
      print "WARNING: variable %s not found in plotfile %s, skipping..." % (variable, dname2)
      continue

    aerror = 0.0

    for level in range(1, plt1.flevel+1):
      for box in range(1, plt1.nboxes(level)+1):

        if refratio == 1:
          err = fdcompare_boxes(
              plt1, comp1, plt2, comp2, level, box, norm)
        else:
          err = fdcompare_boxes_refined(
              plt1, comp1, plt2, comp2, level, box, refratio, norm, diff)

        if norm:
          aerror += err
        else:
          aerror = max(aerror, err)

    if norm:
      errors[variable] = aerror**(1.0/norm)
    else:
      errors[variable] = aerror

  return errors
