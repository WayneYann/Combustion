
import pickle

from collections import namedtuple, defaultdict
from operator import attrgetter
from pprint import pprint
from pylab import *

key = namedtuple('key', [ 'nx', 'cfl', 'nodes' ])

pens = defaultdict(
  lambda: { 'color': 'k', 'marker': 'o' },
  { 0.1: { 'color': 'b', 'marker': 's' },
    0.2: { 'color': 'b', 'marker': 'd' },
    0.3: { 'color': 'b', 'marker': '^' },
    0.4: { 'color': 'b', 'marker': 'v' },
    0.5: { 'color': 'b', 'marker': 'o' },
    0.6: { 'color': 'b', 'marker': 'x' },
    0.7: { 'color': 'b', 'marker': '+' },
    0.8: { 'color': 'b', 'marker': '<' },
    0.9: { 'color': 'b', 'marker': '>' },
    1.0: { 'color': 'b', 'marker': '.' },
    32:  { 'color': 'b', 'marker': 'o' },
    64:  { 'color': 'r', 'marker': '^' },
    128: { 'color': 'k', 'marker': 's' },
    } )


def rekey(d):
  dk = {}
  for idx in d:
    dk[key(idx[0], idx[1], idx[2])] = d[idx]
  return dk


def flameball_stconv_plot(stconv, cflconv):

  # print 'stconv'
  # print '------'

  # figure()

  # nodes = sorted(set([ run.nodes for run in stconv]))
  # cfls  = sorted(set([ run.cfl for run in stconv ]))

  # for j, node in enumerate(nodes):
  #   for cfl in cfls:
  #     runs = [ run for run in stconv 
  #              if run.cfl == cfl and run.nodes == node ]
  #     runs = sorted(runs, key=attrgetter('nx'))

  #     x = np.asarray([ run.nx for run in runs ])
  #     y = np.asarray([ stconv[run] for run in runs ])

  #     print 'nodes: %d, cfl: %f' % (node, cfl)
  #     print '  x:', x
  #     print '  y:', y
  #     print ''
      
  #     subplot(1, 2, j+1)
  #     title('%d nodes' % node)
  #     loglog(x, y, label='cfl %.2f' % cfl, **pens[cfl])
  #     xlabel('nx')
  #     ylabel('l2 error')

  # legend(loc='best')

  print 'cflconv'
  print '-------'

  figure()

  nxs   = sorted(set([ run.nx for run in cflconv]))
  nodes = sorted(set([ run.nodes for run in cflconv]))
  cfls  = sorted(set([ run.cfl for run in cflconv ]))

  for nx in nxs:
    for j, node in enumerate(nodes):
      runs = [ run for run in cflconv 
               if run.nodes == node and run.nx == nx ]
      runs = sorted(runs, key=attrgetter('cfl'))

      x = np.asarray([ run.cfl for run in runs ])
      y = np.asarray([ cflconv[run] for run in runs ])

      print 'nodes: %d, nx: %f' % (node, nx)
      print '  x:', x
      print '  y:', y
      print '  p:', np.log(y[1:] / y[:-1]) / np.log(x[1:] / x[:-1])
      print ''

      subplot(1, len(nodes), j+1)
      title('%d nodes' % node)
      loglog(x, y, label='nx %d' % int(nx), **pens[nx])
      xlabel('dt')
      ylabel('l2 error')

      x = [ x[0], x[-1] ]
      y = [ y[0]*2, y[0]*2*(x[1]/x[0])**(2*node-2) ]
      loglog(x, y, '--k')
      

  legend(loc='best')

if __name__ == '__main__':

  with open('stconv.pkl', 'r') as f:
    errors = pickle.load(f)

  print 'errors from stconv.pkl'
  print '----------------------'
  pprint(errors, indent=2)
  print

  flameball_stconv_plot(rekey(errors['stconv']), rekey(errors['cflconv']))

  show()
  
