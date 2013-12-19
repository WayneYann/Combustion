/* DGNUPLOT related debugging routines */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

#include <MultiFab.H>

#define GPCHKERR(err, msg) if ((err) == -1) { perror(msg); return -1; }
#define MAXWINDOWS 12

FILE *to[MAXWINDOWS] = { 0 };
FILE *fr[MAXWINDOWS] = { 0 };

/*
 * Fork a new gnuplot process and connect pipes.
 */
int dgp_connect(FILE **fto, FILE **ffr)
{
  int err, pid, to[2], fr[2];

  err = pipe(to); GPCHKERR(err, "Unable to create pipe.");
  err = pipe(fr); GPCHKERR(err, "Unable to create pipe.");
  pid = fork();   GPCHKERR(pid, "Unable to fork gnuplot.");

  if (pid) {
    close(to[0]); close(fr[1]);
  } else {
    close(to[1]); close(fr[0]);
    dup2(to[0], 0); close(to[0]);
    dup2(fr[1], 1); close(fr[1]);
    execl("/usr/bin/gnuplot", "gnuplot", NULL);
  }

  *fto = fdopen(to[1], "w");
  *ffr = fdopen(fr[0], "r");

  return 0;
}

/*
 * Send a MultiFab to a gnuplot window.
 */
void dgp_send_mf(MultiFab& U, int level, int comp, int wait)
{
  if (level >= MAXWINDOWS) return;

  if (to[level] == 0) {
    dgp_connect(&to[level], &fr[level]);
  }

  const Box& mbx = U.boxArray().minimalBox();
  fprintf(to[level], "set term wxt %d\n", level);
  fprintf(to[level], "set dgrid3d %d, %d\n", mbx.length(0)+2*U.nGrow(), mbx.length(1)+2*U.nGrow());
  fprintf(to[level], "set pm3d\n");
  fprintf(to[level], "set hidden3d\n");
  fprintf(to[level], "splot '-'\n");
  for (MFIter mfi(U); mfi.isValid(); ++mfi) {
    const Box& bx = U[mfi].box();
    IntVect    sm = bx.smallEnd();
    IntVect    bg = bx.bigEnd();
    for (IntVect p=sm; p<=bg; bx.next(p))
      fprintf(to[level], "%d %d %e\n", p[0], p[1], U[mfi](p,comp));
  }
  fprintf(to[level], "e\n");
  if (wait) {
    fprintf(to[level], "set print\n");
    fprintf(to[level], "pause mouse button3 \"===> paused - button 3 in window %d to release\"\n", level);
    fprintf(to[level], "print ''\n");
    fprintf(to[level], "set print '-'\n");
    fprintf(to[level], "print 'hoser!'\n");
  }
  fflush(to[level]);

  if (wait) {
    char buf[8];
    char *s = fgets(buf, 8, fr[level]);
  }
}
