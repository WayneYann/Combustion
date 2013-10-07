/* DZMQ related debugging routines */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <zmq.h>

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
