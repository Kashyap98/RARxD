
from PySide2.QtWidgets import QWidget, QLineEdit, QVBoxLayout, QComboBox, QLabel, QCheckBox


class CheckBoxWidget(QWidget):
    def __init__(self, message, parent=None):
        super(CheckBoxWidget, self).__init__(parent)
        self.layout = QVBoxLayout()

        self.checkbox = QCheckBox()
        self.checkbox.setChecked(True)

        self.label = QLabel()
        self.label.setText(message)

        self.layout.addWidget(self.label, 1)
        self.layout.addWidget(self.checkbox, 1)

        self.setLayout(self.layout)

    def get_status(self):
        return self.checkbox.isChecked()


class LineEditWidget(QWidget):
    def __init__(self, message, parent=None):
        super(LineEditWidget, self).__init__(parent)
        self.layout = QVBoxLayout()

        self.editbox = QLineEdit()

        self.label = QLabel()
        self.label.setText(message)

        self.layout.addWidget(self.label, 1)
        self.layout.addWidget(self.editbox, 1)

        self.setLayout(self.layout)

    def get_text(self):
        return self.editbox.text()


class ComboBoxWidget(QWidget):
    pass

    def __init__(self, items, message, parent=None):
        super(ComboBoxWidget, self).__init__(parent)
        self.layout = QVBoxLayout()

        self.combobox = QComboBox()
        self.combobox.addItems(items)

        self.label = QLabel()
        self.label.setText(message)

        self.layout.addWidget(self.label, 1)
        self.layout.addWidget(self.combobox, 1)

        self.setLayout(self.layout)

    def get_selection(self):
        return self.combobox.currentText()
