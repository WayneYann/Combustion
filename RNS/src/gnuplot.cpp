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
  fprintf(to[level], "set pm3d corners2color c1 map\n");
  for (MFIter mfi(U); mfi.isValid(); ++mfi) {
    const Box& bx = U[mfi].box();
    IntVect    sm = bx.smallEnd();
    IntVect    bg = bx.bigEnd();

    const int ny = bg[1]-sm[1]+1;
    const int nx = bg[0]-sm[0]+1;

    fprintf(to[level], "splot '-' binary array=%dx%d with pm3d\n", ny, nx);

    IntVect p;
    float r[ny];
    for (int i=0; i<nx; i++) {
      for (int j=0; j<ny; j++) {
	p[0] = sm[0] + i; p[1] = sm[1] + j;
      	r[j] = U[mfi](p,comp);
      }
      fwrite(r, sizeof(float), ny, to[level]);
    }
  }
  if (wait) {
    fprintf(to[level], "set print\n");
    fprintf(to[level], "pause mouse button2 \"===> paused - button 2 in window %d to release\"\n", level);
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
