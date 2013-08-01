
import pickle

from collections import namedtuple, defaultdict
from operator import attrgetter
from pprint import pprint
from pylab import *

class Container():
    pass

pens = { 
  (3, 9):  { 'color': 'r', 'marker': 'o', 'ms': 16, 'mew': 2, 'ls': 'None' },
  (3, 13): { 'color': 'b', 'marker': '^', 'ms': 14, 'mew': 2, 'ls': 'None' },
  }

labels = {
  (3, 9):  'GL 3.9',
  (3, 13): 'GL 3.13',
  }

order = {
  (3, 9):  4,
  (3, 13): 4,
  }


def flameball_mrconv_plot(errors):

  nxs   = sorted(set([ run.nx for run in errors]))
  nodes = sorted(set([ run.nodes for run in errors]))
  dts   = sorted(set([ run.dt for run in errors ]))

  for nx in nxs:
    for j, node in enumerate(nodes):
      runs = [ run for run in errors 
               if run.nodes == node and run.nx == nx ]
      runs = sorted(runs, key=attrgetter('dt'))

      x = np.asarray([ run.dt for run in runs ], np.float64)
      y = np.asarray([ run.error for run in runs ])

      print 'nodes: %s, nx: %f' % (node, nx)
      print '  x:', x
      print '  y:', y
      print '  p:', np.log(y[1:] / y[:-1]) / np.log(x[1:] / x[:-1])
      print ''

      loglog(x, y, label=labels[node], **pens[node])
      
  x = [ x[-1], x[0] ]
  y = [ y[-1], y[-1]*(x[1]/x[0])**order[node] ]
  loglog(x, y, '-k', zorder=-10, label='4th order')

  legend(loc='upper left')
  ylabel('$L_2$ error')
  xlabel('$\Delta t$')
  xlim([ 2.125e-9, 1.1e-8 ])
  savefig('mrconv.eps')




if __name__ == '__main__':

  with open('mrconv.pkl', 'r') as f:
    errors = pickle.load(f)

  print 'errors from mrconv.pkl'
  print '----------------------'
  pprint(errors, indent=2)
  print

  errors = errors['dtconv']

  entries = []
  for entry in errors:
    c = Container()
    c.nx = entry[0]
    c.dt = entry[1]
    c.nodes   = entry[2]
    c.error   = errors[entry][0]
    c.runtime = errors[entry][1][0]
    c.evals   = errors[entry][1][1]
    entries.append(c)

  flameball_mrconv_plot(entries)

  show()
  
