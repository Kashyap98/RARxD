import glob
import os
import sys
import json

import jsonpickle
from PySide2.QtCore import Qt
from PySide2.QtGui import QIcon
from PySide2.QtWidgets import QApplication, QLabel, QLineEdit, QPushButton, QVBoxLayout, QDialog, QHBoxLayout, QWidget, \
    QTabWidget
from gui.dialog_size_util import DialogSizeUtil
from gui.history_panel import HistoryPanel
from gui.input_panel import InputPanel
from gui.run_panel import RunPanel
from helpers.input_data_handler import GlobalApplicationSettings
from helpers.constants import *


class MainDialog(QDialog):

    def __init__(self, q_app, parent=None):
        super(MainDialog, self).__init__(parent)
        self.app = q_app
        self.setWindowFlags(Qt.WindowCloseButtonHint)
        self.setWindowIcon(QIcon(os.path.join(os.getcwd(), "..", "logo.png")))
        self.dialog_size_util = DialogSizeUtil(self.app)
        self.setWindowTitle("Gene Fusion Tool - COMP 383 Project")
        self.resize(self.dialog_size_util.get_window_width(0.67), self.dialog_size_util.get_window_height(0.5))
        self.move(self.dialog_size_util.center_dialog_width(), self.dialog_size_util.center_dialog_height())
        self.history_panel = HistoryPanel(self.app, self)
        self.input_panel = InputPanel(self.app)
        self.run_panel = RunPanel(self.app)
        self.main_layout = QHBoxLayout()
        self.left_panel = QVBoxLayout()
        self.right_panel = QVBoxLayout()

        self.right_panel.addWidget(self.history_panel, 1)
        self.left_panel.addWidget(self.input_panel, 25)
        self.left_panel.addWidget(self.run_panel, 75)

        # self.button.clicked.connect(self.greetings)
        self.main_layout.addLayout(self.left_panel, 72)
        self.main_layout.addLayout(self.right_panel, 28)
        self.setLayout(self.main_layout)
        self.show()

    def view_run(self, folder_path):
        self.run_panel.init_tabs()
        self.log_info_file = open(os.path.join(folder_path, f"{os.path.basename(folder_path)}_log.txt"), "r")
        self.settings_file = open(os.path.join(folder_path, "run_settings.json"), "r")
        self.input_data_file = open(os.path.join(folder_path, "input_data.json"), "r")
        self.run_settings = GlobalApplicationSettings(folder_path, is_writing=False)

        self.run_panel.logs_label.setText(self.log_info_file.read())
        self.run_panel.settings_label.setText(self.build_settings_text())

        self.build_transcript_tabs(folder_path)

    def build_transcript_tabs(self, folder_path):
        transcripts = []
        blast_results = {}
        for folder in os.listdir(folder_path):
            if os.path.isdir(os.path.join(folder_path, folder)):
                transcripts.append(folder)
                transcript_blast_results = []
                for blast_result_file in glob.glob(os.path.join(folder_path, folder, "*.json")):
                    transcript_blast_results.append(jsonpickle.decode(open(blast_result_file, "r").read(), keys=True))
                blast_results[folder] = transcript_blast_results
                self.run_panel.transcript_tab(folder, transcript_blast_results)
        return transcripts, blast_results

    def build_settings_text(self):
        settings_data = jsonpickle.decode(self.settings_file.read(), keys=True)
        input_data = jsonpickle.decode(self.input_data_file.read(), keys=True)
        settings_string = ""
        for key, value in settings_data.items():
            settings_string += f"{key}: {value}\n"
        return settings_string


if __name__ == "__main__":
    app = QApplication(sys.argv)
    main_dialog = MainDialog(app)
    app.exec_()
