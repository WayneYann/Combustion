
import pickle

from collections import namedtuple
from operator import attrgetter
from pylab import *


key = namedtuple('key', [ 'nx', 'cfl', 'nodes' ])


def flameball_stconv_plot():

  with open('stconv.pkl', 'r') as f:
    errors = pickle.load(f)

  nodes = set([ run.nodes for run in errors])
  cfls  = set([ run.cfl for run in errors ])

  for node in nodes:
    for cfl in cfls:
      runs = sorted([ run for run in errors if run.cfl == cfl and run.nodes == node ], key=attrgetter('nx'))
      x = [ run.nx for run in runs ]
      y = [ errors[run] for run in runs ]

      loglog(x, y, label='cfl %.2f, nnodes %d' % (cfl, node))

  legend(loc='best')
  show()






if __name__ == '__main__':
  flameball_stconv_plot()
