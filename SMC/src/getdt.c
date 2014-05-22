
#include <sdc.h>

double get_dt_m(sdc_mrex *mrex, int comp, int m)
{
  sdc_nset*  nset   = mrex->nsets[comp];
  sdc_dtype* nodes  = nset->nodes;
  int        nnodes = nset->nnodes;

  if (m+1 >= nnodes) {
    return nodes[1] - nodes[0];
  }
  return nodes[m+1] - nodes[m];
}
