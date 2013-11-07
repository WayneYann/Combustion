"""Visual debugger."""

import sys
import time
import threading

import zmq
import numpy as np

import pygame

class Acquire(threading.Thread):
    """Listen for data and call 'render_scene' when appropriate.

    This thread creates a ZMQ socket (bound to port 31415) and listens
    for data.  When data is received, the GUI thread's 'render_scene'
    method is called.

    A non-blocking ZMQ receive is used so that we can periodically
    check our 'stop' event and exit if appropriate.  The 'stop' event
    is set by the main thread when the user closes the main window.

    """

    def __init__(self):
        super(Acquire, self).__init__()
        self.stop = threading.Event()
        self.viewer = None


    def recv(self, flag=0):

        try:
            return np.fromstring(self.socket.recv(flag), np.float64)
        except zmq.ZMQError:
            return None


    def run(self):

        self.context = zmq.Context()
        self.socket  = self.context.socket(zmq.PULL)
        self.socket.bind("tcp://*:31415")

        while True:
            if self.stop.is_set():
                return

            rho = self.recv(zmq.NOBLOCK)
            if rho is None:
                time.sleep(0.005)
                continue

            self.viewer.state['rho'] = rho
            self.viewer.state['p']   = self.recv()
            self.viewer.state['v']   = self.recv()
            self.viewer.state['divu']   = self.recv()
            self.viewer.state['macv']  = self.recv()
            self.viewer.state['f_rho'] = self.recv()
            self.viewer.state['f_v']   = self.recv()


class MainWindow():
    """Build main window and render scene when called.

    When the user closes the main window the Acquire thread's stop
    event is set.
    """

    def plot(self, q, x0, y0, w, h, color=(255, 255, 255)):
        x = x0 + w * np.arange(len(q)) / len(q)
        y = y0 + h - h * q / 2.1 - 2
        y[ np.logical_not(np.isfinite(y)) ] = 0
        pygame.draw.aalines(self.screen, color, False, zip(x, y))


    def render(self, q):

        state = self.state

        self.lock.acquire()

        w, h = self.size[0], self.size[1]/4
        r, pad = 5, 10

        x, y = pad, pad
        if 'rho' in state:
            self.plot(state['rho'], x, y, w-2*pad, h-2*pad, color=(0,0,255))

        if 'f_rho' in state:
            self.plot(state['f_rho'], x, y, w-2*pad, h-2*pad, color=(255,255,0))

        if 'p' in state:
            self.plot(state['p'], x, y, w-2*pad, h-2*pad, color=(255,0,0))


        x, y = pad, 200+2*pad
        if 'v' in state:
            self.plot(state['v'], x, y, w-2*pad, h-2*pad, color=(255,0,0))

        if 'divu' in state:
            self.plot(state['divu'], x, y, w-2*pad, h-2*pad, color=(255,0,255))

        if 'macv' in state:
            self.plot(state['macv'], x, y, w-2*pad, h-2*pad, color=(0,255,0))

        if 'f_v' in state:
            self.plot(state['f_v'], x, y, w-2*pad, h-2*pad, color=(255,255,0))

        self.lock.release()


    def __init__(self):

        pygame.init()
        pygame.font.init()

        self.size   = (800, 800)
        self.font   = pygame.font.SysFont('dejavusans', self.size[1]/40)
        self.screen = pygame.display.set_mode(self.size)

        self.state  = {}
        self.lock   = threading.Lock()

        self.acquire = Acquire()
        self.acquire.viewer = self
        self.acquire.start()


    def run(self):

        while True:
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    self.acquire.stop.set()
                    self.acquire.join()
                    sys.exit()

            time.sleep(0.005)

            self.screen.fill((0, 0, 0))
            self.render(None)
            pygame.display.flip()


if __name__ == '__main__':
    frame = MainWindow()
    frame.run()


