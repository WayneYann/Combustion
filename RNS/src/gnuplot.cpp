/* DGNUPLOT related debugging routines */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

#include <map>
#include <iostream>
#include <sstream>

#include <MultiFab.H>

using namespace std;

map<int,FILE*> gp, from;
void dgp_connect(int window);

void dgp_send_mf(MultiFab& U, int level, int comp, int wait)
{
  if (gp.count(level) == 0) dgp_connect(level);

  const Box& mbx = U.boxArray().minimalBox();
  fprintf(gp[level], "set term wxt %d\n", level);
  fprintf(gp[level], "set dgrid3d %d, %d\n", mbx.length(0)+2*U.nGrow(), mbx.length(1)+2*U.nGrow());
  fprintf(gp[level], "set pm3d\n");
  fprintf(gp[level], "set hidden3d\n");
  fprintf(gp[level], "splot '-'\n");
  for (MFIter mfi(U); mfi.isValid(); ++mfi) {
    const Box& bx = U[mfi].box();
    IntVect    sm = bx.smallEnd();
    IntVect    bg = bx.bigEnd();
    for (IntVect p=sm; p<=bg; bx.next(p))
      fprintf(gp[level], "%d %d %e\n", p[0], p[1], U[mfi](p,comp));
  }
  fprintf(gp[level], "e\n");
  if (wait) {
    fprintf(gp[level], "set print\n");
    fprintf(gp[level], "pause mouse button3 \"===> paused - button 3 in window %d to release\"\n", level);
    fprintf(gp[level], "print ''\n");
    fprintf(gp[level], "set print '-'\n");
    fprintf(gp[level], "print 'hoser!'\n");
  }
  fflush(gp[level]);

  if (wait) {
    char buf[8];
    char *s = fgets(buf, 8, from[level]);
  }
}

void dgp_connect(int window)
{
  int err, pid, fdto[2], fdfrom[2];

  err = pipe(fdto);
  if (err == -1) {
    perror("Unable to create pipe");
    return;
  }

  err = pipe(fdfrom);
  if (err == -1) {
    perror("Unable to create pipe");
    return;
  }

  pid = fork();
  if (pid == -1) {
    perror("Unable to fork gnuplot");
    return;
  }

  if (pid) {
    close(fdto[0]);
    close(fdfrom[1]);
  } else {
    close(fdto[1]); close(fdfrom[0]);
    dup2(fdto[0], 0); close(fdto[0]);
    dup2(fdfrom[1], 1); close(fdfrom[1]);
    execl("/usr/bin/gnuplot", "gnuplot", NULL);
  }

  gp[window] = fdopen(fdto[1], "w");
  from[window] = fdopen(fdfrom[0], "r");
}
