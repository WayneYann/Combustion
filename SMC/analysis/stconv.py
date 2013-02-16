
from itertools import product

dt0 = 1e-8
stop_time = dt0 * 200

# list of (nx, dt, nnodes) runs
runs = [ ]

runs.extend(product(
        [ 32 ], 
        [ dt0/4, dt0/2, dt0, 2*dt0 ], 
        [ 3 ]))

# runs.extend(product(
#         [ 128 ],
#         [ dt0/32, dt0/16, dt0/8, dt0/4 ],
#         [ 3, 5 ]))

test_runs = product(
    [ 32, 64 ], [ dt0/8, dt0/4, dt0/2, dt0 ], [ 3, 5 ])

 
min_dt = {
    32: dt0/4,
    64: dt0/4,
    128: dt0/32,
    }
