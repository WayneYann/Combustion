/* DGNUPLOT related debugging routines */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <iostream>
#include <sstream>

#include <MultiFab.H>

FILE *dgp_connect();
FILE* gp = 0;

void dgp_send_mf(MultiFab& U, int level, int comp, int wait)
{
  // if (!zmqctx) zmqctx = dzmq_connect();
  if (!gp) gp = dgp_connect();

  const Box& mbx = U.boxArray().minimalBox();
  fprintf(gp, "set term wxt %d\n", level);
  fprintf(gp, "set dgrid3d %d, %d\n", mbx.length(0)+2*U.nGrow(), mbx.length(1)+2*U.nGrow());
  fprintf(gp, "set pm3d\n");
  fprintf(gp, "set hidden3d\n");
  fprintf(gp, "splot '-'\n");
  for (MFIter mfi(U); mfi.isValid(); ++mfi) {
    const Box& bx = U[mfi].box();
    IntVect    sm = bx.smallEnd();
    IntVect    bg = bx.bigEnd();
    for (IntVect p=sm; p<=bg; bx.next(p))
      fprintf(gp, "%d %d %e\n", p[0], p[1], U[mfi](p,comp));
  }
  fprintf(gp, "e\n");
  if (wait)
    fprintf(gp, "pause mouse button3 \"===> paused - button 3 in window %d to release\"\nprint ''\n", level);
  fflush(gp);

  if (wait) {
    char buf[8];
    printf("===> paused\n");
    fgets(buf, 8, stdin);
  }
}

FILE *dgp_connect()
{
  // if (fork()) {
  //   sleep(1);
    gp = fopen("/tmp/gp", "a");
  // } else {
  //   execl("/usr/bin/gnuplot", "/tmp/gp");
  // }
  return gp;
}
