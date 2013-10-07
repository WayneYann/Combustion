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

tup = "\(" + "([0-9,]+)" + "\)\s*"
parse_header = re.compile("\(" + 3*tup + "\) (\d+)")

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
            # s = np.fromstring(self.socket.recv(flag), np.int32)
            # q = np.fromstring(self.socket.recv())
            # return q.reshape(s, order='F')

            s = self.socket.recv(flag)
            i = s.find('\n')
            hdr = s[:i]
            dat = np.fromstring(s[i+1:], np.float32)

            m = parse_header.search(hdr)
            if m is None:
                return None, None, None
            lo = np.asarray(map(int, m.group(1).split(',')))
            hi = np.asarray(map(int, m.group(2).split(',')))
            dat = dat.reshape(hi-lo+1, order='F')
            return lo, hi, dat

        except zmq.ZMQError:
            return None, None, None


    def run(self):

        self.context = zmq.Context()
        self.socket  = self.context.socket(zmq.PULL)
        self.socket.bind("tcp://*:31415")

        while True:
            if self.stop.is_set():
                return

            lo, hi, u = self.recv(zmq.NOBLOCK)

            if u is None:
                time.sleep(0.005)
                continue

            # v = self.recv()
            # s = self.recv()

            print "u: shape: ", u.shape, ", max: ", u.max(), ", min: ", u.min()
            # print "v: shape: ", v.shape, ", max: ", v.max(), ", min: ", v.min()
            # print "s: shape: ", s.shape, ", max: ", s.max(), ", min: ", s.min()
            # print ""

            if self.viewer is not None:
                self.viewer.render_scene(lo, hi, u)


class MainWindow(QMainWindow):
    """Build main window and render scene when called.

    When the user closes the main window the Acquire thread's stop
    event is set.
    """


    def render_scene(self, lo, hi, u):

        x = np.arange(lo[0], hi[0]+1, 1.0)
        y = np.arange(lo[1], hi[1]+1, 1.0)
        self.ax[0].contourf(x, y, u)

        # self.ax[0].clear()
        # self.ax[0].imshow(u)

        # self.ax[1].clear()
        # self.ax[1].imshow(v) 

        # self.ax[2].clear()
        # self.ax[2].imshow(s) 

        # for ax in self.ax.itervalues():
        #     grid(ax, u.shape[0])

        self.canvas.draw()


    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)

        self.fig = Figure()
        self.ax = {}
        self.ax[0] = self.fig.add_subplot(111)
        
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


