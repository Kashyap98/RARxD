import jsonpickle
from PySide2.QtWidgets import QHBoxLayout, QPushButton, QProgressBar
from PySide2.QtCore import QThread
from PySide2 import QtCore
from gui.input_fields import *
from helpers.setup_helper import handle_setup, GlobalApplicationSettings, os
from helpers.constants import *
from helpers.worker_thread_helper import Worker

PROMOTERS_DICT = {"HSP70": "AF158020.1"}


class InputPanel(QWidget):

    def __init__(self, q_app, history_panel, parent=None):
        super(InputPanel, self).__init__(parent)
        self.app = q_app
        self.layout = QVBoxLayout()
        self.history_panel = history_panel
        self.input_layout = QHBoxLayout()
        self.filter_layout = QHBoxLayout()
        self.submit_layout = QHBoxLayout()
        self.promoters = ["HSP70"]
        self.blast_options = ["LOCAL", "REMOTE"]

        self.run_name_field = LineEditWidget("Run Name:")
        self.accession_field = LineEditWidget("Ensembl Accession ENSDARGXXXXXXXXXXX:")
        self.promoter_combobox = ComboBoxWidget(self.promoters, "Promoter:")
        self.blast_combobox = ComboBoxWidget(self.blast_options, "BLAST Option:")

        self.filter_complete = CheckBoxWidget("Filter Complete Transcripts?")
        self.filter_protein_coding = CheckBoxWidget("Filter Protein Coding Transcripts?")
        self.filter_start_codon = CheckBoxWidget("Filter Start Codon Missing Transcripts?")
        self.filter_stop_codon = CheckBoxWidget("Filter Stop Codon Missing Transcripts?")

        self.start_button = QPushButton()
        self.progress_bar = QProgressBar()
        self.start_button.setText("Start Run")
        self.start_button.clicked.connect(self.start_run)

        self.input_layout.addWidget(self.run_name_field)
        self.input_layout.addWidget(self.accession_field)
        self.input_layout.addWidget(self.promoter_combobox)
        self.input_layout.addWidget(self.blast_combobox)

        self.filter_layout.addWidget(self.filter_complete)
        self.filter_layout.addWidget(self.filter_protein_coding)
        self.filter_layout.addWidget(self.filter_start_codon)
        self.filter_layout.addWidget(self.filter_stop_codon)

        self.submit_layout.addWidget(self.progress_bar)
        self.submit_layout.addWidget(self.start_button)

        self.layout.addLayout(self.input_layout)
        self.layout.addLayout(self.filter_layout)
        self.layout.addLayout(self.submit_layout)
        self.setLayout(self.layout)

        self.worker_thread = QThread()

    def start_run(self):
        run_name = self.run_name_field.get_text().replace(" ", "")
        ensembl_id = self.accession_field.get_text().replace(" ", "")

        if run_name != "" and ensembl_id != "":
            promoter_id = self.promoter_combobox.get_selection()
            blast_type = self.blast_combobox.get_selection()

            filter_complete = self.filter_complete.get_status()
            filter_protein_coding = self.filter_protein_coding.get_status()
            filter_start_codon = self.filter_start_codon.get_status()
            filter_stop_codon = self.filter_stop_codon.get_status()

            global_data, global_settings, folder_path = handle_setup(gui=True, folder_path=run_name)

            output_settings = GlobalApplicationSettings(folder_path)
            output_settings.add_to_settings(FILTER_TRANSCRIPT_INCOMPLETE, filter_complete)
            output_settings.add_to_settings(FILTER_TRANSCRIPT_PROTEIN_CODING, filter_protein_coding)
            output_settings.add_to_settings(FILTER_TRANSCRIPT_START_CODON, filter_start_codon)
            output_settings.add_to_settings(FILTER_TRANSCRIPT_STOP_CODON, filter_stop_codon)
            output_settings.export_as_json()

            output_data = {TARGET_SEQUENCE: ensembl_id, PROMOTER: PROMOTERS_DICT[promoter_id], BLAST_TYPE: blast_type}
            self.export_input_data(folder_path, output_data)

            controller_path = os.path.join(os.getcwd(), "controller.py")

            worker = Worker(controller_path, run_name)
            worker.moveToThread(self.worker_thread)
            self.worker_thread.started.connect(worker.execute_run)
            self.worker_thread.finished.connect(self.refresh_data)
            self.worker_thread.start()

    @QtCore.Slot()
    def refresh_data(self):
        print("finished")
        self.history_panel.set_data()

    def export_input_data(self, folder_path, output_data):
        with open(os.path.join(folder_path, INPUT_DATA_FILE_NAME), "w") as output_file:
            output_file.write(jsonpickle.encode(output_data, keys=True))

