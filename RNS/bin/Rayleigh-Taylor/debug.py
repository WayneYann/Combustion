"""Visual debugger."""

import sys
import time
import threading

import zmq
import numpy as np

from PySide.QtCore import *
from PySide.QtGui import *

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import re

tup = "\(" + "([-0-9,]+)" + "\)\s*"
parse_header = re.compile("\(" + 3*tup + "\) (\d+)")

ng = 4

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
            s = self.socket.recv(flag)
            level = s[0]
            s = s[1:]
            i = s.find('\n')
            hdr = s[:i]
            dat = np.fromstring(s[i+1:])

            m = parse_header.search(hdr)
            if m is None:
                print "error parsing header: ", hdr
                return None, None, None, None
            lo = np.asarray(map(int, m.group(1).split(',')))
            hi = np.asarray(map(int, m.group(2).split(',')))
            dat = dat.reshape(hi-lo+1, order='F')
            return int(level), lo, hi, dat

        except zmq.ZMQError:
            return None, None, None, None


    def run(self):

        self.context = zmq.Context()
        self.socket  = self.context.socket(zmq.PULL)
        self.socket.bind("tcp://*:31415")

        while True:
            if self.stop.is_set():
                return

            level, lo, hi, u = self.recv(zmq.NOBLOCK)

            if u is None:
                time.sleep(0.005)
                continue

            # v = self.recv()
            # s = self.recv()

            print "u: shape: ", u.shape, \
              ", max: ", u[ng:-ng,ng:-ng].max(), \
              ", min: ", u[ng:-ng,ng:-ng].min(), \
              ", lo/hi:", lo, hi
            # print "v: shape: ", v.shape, ", max: ", v.max(), ", min: ", v.min()
            # print "s: shape: ", s.shape, ", max: ", s.max(), ", min: ", s.min()
            # print ""

            if self.viewer is not None:
                self.viewer.render_scene(level, lo, hi, u)


class MainWindow(QMainWindow):
    """Build main window and render scene when called.

    When the user closes the main window the Acquire thread's stop
    event is set.
    """


    def render_scene(self, level, lo, hi, u):

        x = np.arange(lo[0], hi[0]+1, 1.0)
        y = np.arange(lo[1], hi[1]+1, 1.0)
        if level > 0:
            x = x / 2**level
            y = y / 2**level
        self.ax[level].contourf(x, y, u.transpose(), vmin=1.0, vmax=2.0, antialias=False)

        # self.ax[0].cla()
        # self.ax[0].matshow(u.transpose())


        self.canvas.draw()


    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)

        self.fig = Figure()
        self.ax = {}
        self.ax[0] = self.fig.add_subplot(211)
        self.ax[1] = self.fig.add_subplot(212)
        
        self.canvas = FigureCanvas(self.fig)
        self.setCentralWidget(self.canvas)

        self.acquire = Acquire()
        self.acquire.viewer = self
        self.acquire.start()


    def closeEvent(self, event):

        self.acquire.stop.set()
        self.acquire.join()
        event.accept()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    frame = MainWindow()
    frame.show()
    sys.exit(app.exec_())


