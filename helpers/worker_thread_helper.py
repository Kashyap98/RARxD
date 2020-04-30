import os
import threading
import time

from PySide2 import QtCore, QtWidgets
from PySide2.QtCore import Slot
from PySide2.QtWidgets import QProgressDialog, QApplication


class RunHelperThread(threading.Thread):

    def __init__(self, command=""):
        threading.Thread.__init__(self)
        self.command = command
        self.setDaemon(True)
        self.start()

    def run(self):
        os.system(self.command)
