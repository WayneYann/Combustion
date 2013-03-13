
import pickle

from collections import namedtuple, defaultdict
from operator import attrgetter
from pprint import pprint
from pylab import *

key = namedtuple('key', [ 'nx', 'dt', 'nodes' ])

pens = defaultdict(
  lambda: { 'color': 'k', 'marker': 'o' },
  { 'gl3': { 'color': 'b', 'marker': 's' },
    'rk3': { 'color': 'k', 'marker': 'd' },
    'gl3-5': { 'color': 'r', 'marker': 'd' },
    'gl3-9': { 'color': 'r', 'marker': 's' },
    } )


def rekey(d):
  dk = {}
  for idx in d:
    dk[key(idx[0], idx[1], idx[2])] = d[idx]
  return dk


order = {
  'rk3': 3,
  'gl3': 4,
  'gl3-5': 4,
  'gl3-9': 4,
  }


def flameball_mrconv_plot(dtconv):

  # figure(1)

  nxs   = sorted(set([ run.nx for run in dtconv]))
  nodes = sorted(set([ run.nodes for run in dtconv]))
  dts   = sorted(set([ run.dt for run in dtconv ]))

  for nx in nxs:
    for j, node in enumerate(nodes):
      runs = [ run for run in dtconv 
               if run.nodes == node and run.nx == nx ]
      runs = sorted(runs, key=attrgetter('dt'))

      x = np.asarray([ dtconv[run][1][1] for run in runs ])
      y = np.asarray([ dtconv[run][0] for run in runs ])

      # x = np.asarray([ run.dt for run in runs ])
      # y = np.asarray([ dtconv[run][0] for run in runs ])

      print 'nodes: %s, nx: %f' % (node, nx)
      print '  x:', x
      print '  y:', y
      print '  p:', np.log(y[1:] / y[:-1]) / np.log(x[1:] / x[:-1])
      print ''

      # subplot(1, len(nodes), j+1)
      title(node)
      loglog(x, y, label=node, **pens[node])
      xlabel('dt')
      ylabel('l2 error')

      # x = [ x[0], x[-1] ]
      # y = [ y[0]*2, y[0]*2*(x[1]/x[0])**order[node] ]
      # loglog(x, y, '--k')

  legend(loc='best')



if __name__ == '__main__':

  with open('mrconv.pkl', 'r') as f:
    errors = pickle.load(f)

  print 'errors from mrconv.pkl'
  print '----------------------'
  pprint(errors, indent=2)
  print

  flameball_mrconv_plot(rekey(errors['dtconv']))

  show()
  
