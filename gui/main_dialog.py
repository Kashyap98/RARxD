import sys

from PySide2.QtCore import Qt
from PySide2.QtWidgets import QApplication, QLabel, QLineEdit, QPushButton, QVBoxLayout, QDialog, QHBoxLayout, QWidget, \
    QTabWidget
from gui.dialog_size_util import DialogSizeUtil
from gui.history_panel import HistoryPanel
from gui.input_panel import InputPanel


class RunPanel(QWidget):

    def __init__(self, q_app, parent=None):
        super(RunPanel, self).__init__(parent)
        self.app = q_app
        layout = QVBoxLayout()
        self.edit = QLineEdit("Write my name here..")
        self.button = QPushButton("Show Greetings")
        layout.addWidget(self.edit)
        layout.addWidget(self.button)
        # self.tab_widget = QTabWidget()
        # self.tab_widget.addTab(self.history_panel, "History")
        self.setLayout(layout)


class MainDialog(QDialog):

    def __init__(self, q_app, parent=None):
        super(MainDialog, self).__init__(parent)
        self.app = q_app
        self.setWindowFlags(Qt.WindowCloseButtonHint)
        self.dialog_size_util = DialogSizeUtil(self.app)
        self.setWindowTitle("Gene Fusion Tool - COMP 383 Project")
        self.resize(self.dialog_size_util.get_window_width(0.67), self.dialog_size_util.get_window_height(0.5))
        self.move(self.dialog_size_util.center_dialog_width(), self.dialog_size_util.center_dialog_height())
        self.history_panel = HistoryPanel(self.app)
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


if __name__ == "__main__":
    app = QApplication(sys.argv)
    main_dialog = MainDialog(app)
    app.exec_()
