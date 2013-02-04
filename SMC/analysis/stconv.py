
from itertools import product

dt0 = 1e-7
stop_time = dt0 * 10

# list of (nx, dt, nnodes) runs
runs = [ ]

# runs.extend(product(
#         [ 32, 64 ], 
#         [ dt0/8, dt0/4, dt0/2, dt0, 2*dt0, 4*dt0 ], 
#         [ 3, 5 ]))

runs.extend(product(
        [ 128 ],
        [ dt0/32, dt0/16, dt0/8, dt0/4 ],
        [ 3, 5 ]))
