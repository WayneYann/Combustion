/* DZMQ related debugging routines */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>

#include <MultiFab.H>
#include <zmq.h>

void *zmqctx = 0;
void *dzmq_connect();
void dzmq_send_buf(void *ptr, const char *buf, int n);

void dzmq_send_mf(MultiFab& U, int level, int comp, int wait)
{
  if (!zmqctx) zmqctx = dzmq_connect();

  int n = 0;
  for (MFIter mfi(U); mfi.isValid(); ++mfi) { n++; }

  int i = 0;
  for (MFIter mfi(U); mfi.isValid(); ++mfi) {
    std::ostringstream buf;
    buf << level << ", " << i++ << ", " << n << ";";
    U[mfi].writeOn(buf, comp, 1);
    dzmq_send_buf(zmqctx, buf.str().c_str(), buf.str().length());
  }

  if (wait) {
    char buf[8];
    printf("===> paused\n");
    fgets(buf, 8, stdin);
  }
}

typedef struct {
  void *context;
  void *socket;
} zmq_ctx_t;


void dzmq_free(void *data, void *hint)
{
  free(data);
}


void *dzmq_connect()
{
  zmq_ctx_t *ctx;
  int rc;

  ctx = (zmq_ctx_t *) malloc(sizeof(zmq_ctx_t));
  if (ctx == NULL) {
    perror("Error: zmq allocate:");
    return NULL;
  }

  ctx->context = zmq_init(1);
  if (ctx->context == NULL) {
    fprintf(stderr, "Error: zmq create: %s\n",
	    zmq_strerror(errno));
  }

  ctx->socket = zmq_socket(ctx->context, ZMQ_PUSH);
  rc = zmq_connect(ctx->socket, "tcp://127.0.0.1:31415");

  if (rc != 0) {
    fprintf(stderr, "Error: zmq connect: %s\n",
	    zmq_strerror(errno));
    // XXX: tear down properly
    return NULL;
  }

  return ctx;
}

void dzmq_close(void *ptr)
{
  zmq_ctx_t *ctx;
  ctx = (zmq_ctx_t *) ptr;

  zmq_close(ctx->socket);
  zmq_term(ctx->context);
}


void dzmq_send_size(void *ptr, int nx, int ny)
{
  zmq_ctx_t *ctx;
  zmq_msg_t msg;
  int *qp;

  ctx = (zmq_ctx_t *) ptr;

  qp = (int *) malloc(sizeof(int)*2);
  qp[0] = nx;
  qp[1] = ny;

  zmq_msg_init_data(&msg, qp, sizeof(int)*2, dzmq_free, NULL);
  zmq_send(ctx->socket, &msg, 0);

  /* zmq_sendmsg(ctx->socket, &msg, 0); */
}

void dzmq_send(void *ptr, double *q, int n)
{
  zmq_ctx_t *ctx;
  zmq_msg_t msg;
  double *qp;

  ctx = (zmq_ctx_t *) ptr;

  qp = (double *) malloc(sizeof(double)*n);
  memcpy(qp, q, sizeof(double)*n);

  zmq_msg_init_data(&msg, qp, sizeof(double)*n, dzmq_free, NULL);
  zmq_send(ctx->socket, &msg, 0);

  /* zmq_sendmsg(ctx->socket, &msg, 0); */
}

void dzmq_send_buf(void *ptr, const char *buf, int n)
{
  zmq_ctx_t *ctx;
  zmq_msg_t msg;
  char *cp;

  ctx = (zmq_ctx_t *) ptr;

  cp = (char *) malloc(n);
  memcpy(cp, buf, n);

  zmq_msg_init_data(&msg, cp, n, dzmq_free, NULL);
  zmq_send(ctx->socket, &msg, 0);

  /* zmq_sendmsg(ctx->socket, &msg, 0); */
}
