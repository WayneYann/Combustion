from math import *
import numpy as np
from scipy.interpolate import interp1d, griddata
from scipy.integrate import cumtrapz, quad
import pdb

import pde
import ml

FT = 1.0
Nx = 10000

a = -1.0
eps = 20.0
r = -50.0

soln = pde.solve_pde(a, eps, r, Nx=Nx, FT=FT, plot_it=True)

np.savetxt('fine', soln.T)
