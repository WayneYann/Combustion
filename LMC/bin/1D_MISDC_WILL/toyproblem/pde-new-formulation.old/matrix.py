from __future__ import division
import numpy as np
import sys
import pdb
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix
from scipy.interpolate import interp1d
from math import *

np.set_printoptions(linewidth=10000)

Nx = 15
h = 1.0/Nx

# identity
I = np.eye(Nx-2)
# create the differentiation matrices
#Dx = (np.diag([1]*(Nx-3),1) - np.diag([1]*(Nx-3),-1))/(2*h)
Dx = (np.diag([1.0/12.0]*(Nx-4),-2)
    + np.diag([-2.0/3.0]*(Nx-3),-1)
    + np.diag([2.0/3.0]*(Nx-3) , 1)
    + np.diag([-1.0/12.0]*(Nx-4),2))/h
Dx[0] = np.array([-5.0/6.0, 3.0/2.0, -1.0/2.0, 1.0/12.0]+[0]*(Nx-6))
Dx[-1] = np.array([0]*(Nx-6) + [-1.0/12.0, 1.0/2.0, -3.0/2.0, 5.0/6.0])

print Dx
