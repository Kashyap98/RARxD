from PySide2.QtWidgets import QHBoxLayout
from gui.input_fields import *


class InputPanel(QWidget):

    def __init__(self, q_app, parent=None):
        super(InputPanel, self).__init__(parent)
        self.app = q_app
        self.layout = QVBoxLayout()
        self.input_layout = QHBoxLayout()
        self.filter_layout = QHBoxLayout()
        self.submit_layout = QHBoxLayout()
        self.promoters = ["HSP70"]

        self.run_name_field = LineEditWidget("Run Name:")
        self.accession_field = LineEditWidget("Ensembl Accession:")
        self.promoter_combobox = ComboBoxWidget(self.promoters, "Promoter:")

        self.filter_complete = CheckBoxWidget("Filter Complete Transcripts?")
        self.filter_protein_coding = CheckBoxWidget("Filter Protein Coding Transcripts?")
        self.filter_start_codon = CheckBoxWidget("Filter Start Codon Missing Transcripts?")
        self.filter_stop_codon = CheckBoxWidget("Filter Stop Codon Missing Transcripts?")

        self.input_layout.addWidget(self.run_name_field)
        self.input_layout.addWidget(self.accession_field)
        self.input_layout.addWidget(self.promoter_combobox)

        self.filter_layout.addWidget(self.filter_complete)
        self.filter_layout.addWidget(self.filter_protein_coding)
        self.filter_layout.addWidget(self.filter_start_codon)
        self.filter_layout.addWidget(self.filter_stop_codon)

        self.layout.addLayout(self.input_layout)
        self.layout.addLayout(self.filter_layout)
        self.layout.addLayout(self.submit_layout)
        self.setLayout(self.layout)
