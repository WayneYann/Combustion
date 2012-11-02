
from __future__ import with_statement

import pickle
from pyboxlib.utils import auto_compare


comps = {}

reference = 'cflconv-sdc/nnodes9_cflfac0.100000'

comps['sdc'] = auto_compare('cflconv-sdc', [ ('sdc_nnodes', 'SDC'), ('cflfac', 'CFL') ], reference)
comps['rk']  = auto_compare('cflconv-rk',  [ ('cflfac', 'CFL') ], reference)

with open('comps.pkl', 'w') as f:
    pickle.dump(comps, f)

