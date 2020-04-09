import glob
import os

from PySide2.QtWidgets import QWidget, QLineEdit, QVBoxLayout, QPushButton, QTableWidgetItem, QListWidget, \
    QListWidgetItem, QTableWidget
from PySide2 import QtCore, QtWidgets
from PySide2.QtCore import Qt, Slot, Signal
from cleaner import delete_folders, remove_folder

from helpers.constants import *


class DeleteRunWidget(QPushButton):

    def __init__(self, folder_path, parent=None):
        super(DeleteRunWidget, self).__init__(parent)
        self.folder_path = folder_path
        self.parent = parent
        self.button = QPushButton()
        self.button.setText("Delete Run")
        self.button.clicked.connect(self.delete_folder)

    def delete_folder(self):
        if remove_folder(self.folder_path):
            self.parent.setData()


class ViewRunWidget(QPushButton):

    def __init__(self, folder_path, parent=None):
        super(ViewRunWidget, self).__init__(parent)
        self.folder_path = folder_path
        self.parent = parent
        self.button = QPushButton()
        self.button.setText("View Run")
        self.button.clicked.connect(self.view_folder)

    def view_folder(self):
        self.parent.parent.view_run(self.folder_path)


def get_folder_names():
    all_folders = glob.glob(os.path.join(os.getcwd(), "..", f"{FOLDER_PREFIX}*"))
    output_dataset = {}
    for folder in all_folders:
        output_dataset[os.path.basename(folder)[4:]] = folder
    return output_dataset


class HistoryPanel(QWidget):

    def __init__(self, q_app, parent=None):
        super(HistoryPanel, self).__init__(parent)
        self.parent = parent
        self.app = q_app
        self.data = get_folder_names()
        self.table = QTableWidget(0, 3, self)
        self.table.setHorizontalHeaderLabels(('Run Name', 'View Run', 'Delete Run'))

        self.layout = QVBoxLayout()
        self.setData()

        self.delete_all_button = QPushButton()
        self.delete_all_button.setText("Delete All Runs")
        self.delete_all_button.clicked.connect(self.delete_all_folders)

        self.layout.addWidget(self.delete_all_button, 1)
        self.layout.addWidget(self.table, 1)
        self.setLayout(self.layout)

    def delete_all_folders(self):
        all_folders = glob.glob(os.path.join(os.getcwd(), "..", f"{FOLDER_PREFIX}*"))
        delete_folders(all_folders)
        self.setData()

    def setData(self):
        count = 0
        self.data = get_folder_names()
        self.table.setRowCount(len(self.data))
        for name, path in self.data.items():
            if name != "":
                name_item = QTableWidgetItem(name)
                self.table.setItem(count, 0, name_item)
                self.table.setCellWidget(count, 1, ViewRunWidget(path, parent=self).button)
                self.table.setCellWidget(count, 2, DeleteRunWidget(path, parent=self).button)
                count += 1


