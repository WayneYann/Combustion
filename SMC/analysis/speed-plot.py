
import pickle
import numpy as np

from pylab import *
from pprint import pprint

class Container():
    pass


def flamebox_speed_plot(speed):

    schemes = set([ x.scheme for x in speed ])

    plotkwargs = {
        'gl3': { 'marker': 'o', 'color': 'black', 'label': 'SR GL 3' },
        'gl3.9': { 'marker': 's', 'color': 'blue', 'label': 'MR GL 3, 9' },
        'gl3.13': { 'marker': '^', 'color': 'red' , 'label': 'MR GL 3, 13' },
        }

    # runtime vs dt
    figure()
    for scheme in schemes:
        x, y = np.asarray(sorted([ (r.dt, r.runtime) for r in speed 
                                   if r.scheme == scheme ])).transpose()
        semilogx(x, y, **plotkwargs[scheme])

    ylim([0, 1500])
    legend()
    xlabel('dt')
    ylabel('runtime')
    savefig('mrsdc_runtime_vs_dt.eps')


    # adv/diff evals vs dt
    figure()
    for scheme in schemes:
        x, y = np.asarray(sorted([ (r.dt, r.evals) for r in speed 
                                   if r.scheme == scheme ])).transpose()
        semilogx(x, y, **plotkwargs[scheme])

    ylim([0, 3000])
    legend()
    xlabel('dt')
    ylabel('number of adv/diff evals')
    savefig('mrsdc_nfevals_vs_dt.eps')

    
    # # runtime vs error
    # figure()
    # for scheme in schemes:
    #     x, y = np.asarray(sorted([ (r.error, r.runtime) for r in speed 
    #                                if r.scheme == scheme ])).transpose()
    #     semilogx(x, y, **plotkwargs[scheme])

    # ylim([0, 3000])
    # legend()
    # xlabel('error')
    # ylabel('runtime')

    # # adv/diff evals vs error
    # figure()
    # for scheme in schemes:
    #     x, y = np.asarray(sorted([ (r.error, r.evals) for r in speed 
    #                                if r.scheme == scheme ])).transpose()
    #     semilogx(x, y, **plotkwargs[scheme])

    # ylim([0, 3000])
    # legend()
    # xlabel('error')
    # ylabel('number of adv/diff evals')

    # # error vs dt
    # print 'error vs dt:'

    # figure()
    # for scheme in schemes:
    #     x, y = np.asarray(sorted([ (r.dt, r.error) for r in speed 
    #                                if r.scheme == scheme ])).transpose()
    #     loglog(x, y, **plotkwargs[scheme])

    #     print scheme, '  p:', np.log(y[1:] / y[:-1]) / np.log(x[1:] / x[:-1])

    # legend()
    # xlabel('dt')
    # ylabel('error')


if __name__ == '__main__':

  with open('speed.pkl', 'r') as f:
    speed = pickle.load(f)

  entries = []
  for entry in speed:
      c = Container()
      c.scheme  = entry[0]
      c.dt      = entry[1]
      c.error   = speed[entry][0]
      c.runtime = speed[entry][1][0]
      c.evals   = speed[entry][1][1]
      entries.append(c)

  print 'contents of speed.pkl'
  print '---------------------'
  pprint(speed, indent=2)
  print

  flamebox_speed_plot(entries)

  show()

