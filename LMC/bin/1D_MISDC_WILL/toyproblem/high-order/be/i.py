import numpy as np
from math import *

# this is unfortunately hard coded according to the number of substeps...
# Gauss-Lobatto quadrature points
quad_pts = [0, 0.276393, 0.723607, 1]
# number of advection substeps
N_A = len(quad_pts)
# quadrature weights computed using Mathematica
weights = [[0.2206011330,   0.3793988670, -0.06781472846, 0.02060113296],
           [-0.07453559925, 0.5217491947,  0.5217491947, -0.07453559925],
           [0.02060113296, -0.06781472846, 0.3793988670,  0.2206011330]]
# perform our Gaussian sub-interval integration
def int_4(f, m, l):
    return np.dot(f, weights[m])*l*0.5


f = [1.,1.69137,2.8117,3.50517]

print int_4(f, 0, .1)
print int_4(f, 1, .1)
print int_4(f, 2, .1)
