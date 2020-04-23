import os
import time

from PySide2 import QtCore, QtWidgets
from PySide2.QtCore import Slot
from PySide2.QtWidgets import QProgressDialog


class Worker(QtCore.QObject):

    def __init__(self, controller_path, run_name, parent=None):
        super(self.__class__, self).__init__(parent)
        self.controller_path = controller_path
        self.run_name = run_name

    @QtCore.Slot()
    def execute_run(self):
        os.system(f"python {self.controller_path} --name {self.run_name}")


class ProgressBarThread(QtCore.QObject):

    def __init__(self, max_value, run_name):
        super().__init__()
        self.max_value = max_value
        self.run_name = run_name
        self.progress_value = 0

        self.progress_dialog = QProgressDialog(f"Performing run: {self.run_name}", "Cancel",
                                               self.progress_value, self.max_value)
        self.progress_dialog.setWindowTitle("Run in progress...")
        self.progress_dialog.setAutoClose(True)

    @QtCore.Slot()
    def open_dialog(self):
        self.progress_dialog.forceShow()

    def increment_progress_bar(self):
        self.progress_value += 1
        self.progress_dialog.setValue(self.progress_value/9)
        QtWidgets.QApplication.processEvents()
