import sys
import threading
import time

from PySide2.QtWidgets import QDialog, QApplication, QLabel, QVBoxLayout, QProgressBar

from gui.dialog_size_util import DialogSizeUtil


class ProgressBarUtil(QDialog):

    def __init__(self, app=QApplication.instance(), max_completed=100, title="Performing Run"):
        super(ProgressBarUtil, self).__init__()
        self.app = app
        self.setWindowTitle(title)
        self.completed = 0
        self.max_completed = max_completed
        self.completion_rate = self.calculate_completion_rate()
        self.resize_dialog()
        self.label = QLabel()
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(self.dialog_width * 0.05, self.dialog_height * 0.05,
                                       self.dialog_width * 0.05, self.dialog_height * 0.05)

        self.progress = QProgressBar(self)
        self.progress.setGeometry(self.dialog_width * 0.05, self.dialog_height * 0.05,
                                  self.dialog_width * 0.90, self.dialog_height * 0.10)

        main_layout.addWidget(self.label)
        main_layout.addWidget(self.progress)

        self.show()

    def resize_dialog(self):
        self.dialog_size_utils = DialogSizeUtil(self.app)
        self.dialog_width = self.dialog_size_utils.get_window_width(0.6)
        self.dialog_height = self.dialog_size_utils.get_window_height(0.3)
        self.resize(self.dialog_width, self.dialog_height * 0.35)
        self.move((self.dialog_size_utils.dialog_width-self.dialog_width)/2,
                  (self.dialog_size_utils.dialog_height-self.dialog_height)/2)

    def calculate_completion_rate(self):
        try:
            completion_rate = round(100/self.max_completed)
        except ZeroDivisionError:
            completion_rate = 1
        return completion_rate

    def increment_progress(self, new_label_text=""):
        self.label.setText(new_label_text)
        if self.max_completed != 0:
            if self.completed < 100:
                self.completed += 1
                self.progress.setValue(self.completed * self.completion_rate)
            self.app.processEvents()
            self.dialog_close_check()

    def dialog_close_check(self):
        if self.completed >= self.max_completed:
            self.close()
            return True
        return False

    def show_progress_bar(self):
        self.keep_processing_events()

    def keep_processing_events(self):
        result = self.dialog_close_check()
        if not result:
            time.sleep(0.25)
            self.keep_processing_events()


class ProgressBarThread(threading.Thread):

    def __init__(self, app=QApplication.instance(), max_completed=100, title="Starting Run"):
        threading.Thread.__init__(self)
        self.app = app
        self.max_completed = max_completed
        self.title = title
        self.setDaemon(True)
        self.progress_bar = ProgressBarUtil(self.app.instance(), self.max_completed, self.title)
        self.start()

    def increment_progress(self, new_label_text=""):
        self.progress_bar.increment_progress(new_label_text)

    def run(self):
        self.progress_bar.show_progress_bar()


# progress bar test code
if __name__ == "__main__":
    app = QApplication(sys.argv)
    GUI = ProgressBarThread(app=app, max_completed=10)
    for i in range(10):
        time.sleep(0.5)
        GUI.increment_progress()
