
#include <SDCAmr.H>
#include "RNS.H"

using namespace std;

sdc_sweeper* rns_sdc_build_level(int lev)
{
  int nnodes0 = 3;
  int trat    = 2;
  int nnodes  = 1 + (nnodes0 - 1) * ((int) pow(trat, lev));

  sdc_nodes* nodes = sdc_nodes_create(nnodes, SDC_GAUSS_LOBATTO);
  sdc_imex*  imex  = sdc_imex_create(nodes, sdc_feval, NULL, NULL);

  sdc_nodes_destroy(nodes);
  sdc_imex_setup(imex, NULL, NULL);

  return (sdc_sweeper*) imex;
}
