import os
import platform
import subprocess

from PySide2.QtCore import Qt
from PySide2.QtWidgets import QVBoxLayout, QWidget, QTabWidget, QLabel, QScrollArea


class BlastResultView(QWidget):

    def __init__(self, result, parent=None):
        super(BlastResultView, self).__init__(parent)
        self.layout = QVBoxLayout()
        self.ncbi_link = QLabel(f"'<a href=\'https://www.ncbi.nlm.nih.gov/nucleotide/{result.accession}\'>'"
                                f"Click this link to view NCBI Page'</a>'")
        self.ncbi_link.setOpenExternalLinks(True)
        self.name_label = QLabel(f"Accession: {result.accession} - {result.record['description']}"
                                 f" - {result.record['organism']}")
        self.scores_label = QLabel(f"Bitscore: {result.bitscore} \n {result.match_length}")

        self.name_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.scores_label.setTextInteractionFlags(Qt.TextSelectableByMouse)

        self.layout.addWidget(self.ncbi_link)
        self.layout.addWidget(self.name_label)
        self.layout.addWidget(self.scores_label)

        self.setLayout(self.layout)


class TranscriptResultView(QWidget):

    def __init__(self, result, parent=None):
        super(TranscriptResultView, self).__init__(parent)
        self.path = os.path.join(result.global_data.folder_path, result.transcript.transcript_name)
        self.layout = QVBoxLayout()
        self.ensembl_link = QLabel(f" '<a href=\'http://useast.ensembl.org/Danio_rerio/Gene/Summary?g="
                                   f"{result.transcript.gene_id}\'>' Click this link to view Ensembl Page'</a>'")
        self.ensembl_link.setOpenExternalLinks(True)
        self.transcript_info = QLabel(f"Gene Name: {result.transcript.gene_name} \n"
                                      f"Transcript Name: {result.transcript.transcript_name} \n"
                                      f"Ensembl ID: {result.transcript.transcript_id} \n"
                                      f"Strand: {result.transcript.strand} \n")
        self.transcript_data = QLabel(f"Biotype: {result.transcript.biotype} \n"
                                      f"Complete: {result.transcript.complete} \n"
                                      f"Contains start codon: {result.transcript.contains_start_codon} \n"
                                      f"Contains stop codon: {result.transcript.contains_stop_codon} \n")
        self.folder_open = QLabel("<a style='color: blue;'> Click this to open run folder </a>")
        self.folder_open.mousePressEvent = self.open_file
        self.transcript_sequences = QLabel(f"Coding Sequence: {result.transcript.coding_sequence} \n"
                                           f"5' UTR: {result.transcript.five_prime_utr_sequence} \n"
                                           f"3' UTR: {result.transcript.three_prime_utr_sequence}")

        self.transcript_info.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.transcript_data.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.transcript_sequences.setTextInteractionFlags(Qt.TextSelectableByMouse)

        self.layout.addWidget(self.ensembl_link)
        self.layout.addWidget(self.transcript_info)
        self.layout.addWidget(self.transcript_data)
        self.layout.addWidget(self.folder_open)
        self.layout.addWidget(self.transcript_sequences)

        self.setLayout(self.layout)

    def open_file(self, event):
        if platform.system() == "Windows":
            subprocess.Popen(["explorer", "/select,", self.path])
        elif platform.system() == "Darwin":
            subprocess.Popen(["open", self.path])
        else:
            subprocess.Popen(["xdg-open", self. path])


class RunPanel(QWidget):

    def __init__(self, q_app, parent=None):
        super(RunPanel, self).__init__(parent)
        self.app = q_app
        self.layout = QVBoxLayout()

        self.tabs_widget = QTabWidget()

        self.layout.addWidget(self.tabs_widget)
        self.setLayout(self.layout)

        self.init_tabs()

    def init_tabs(self):
        while self.tabs_widget.count() > 0:
            self.tabs_widget.removeTab(0)
        self.tabs_widget.addTab(self.logs_tab(), "Run Logs")
        self.tabs_widget.addTab(self.settings_tab(), "Run Settings")

    def logs_tab(self):
        logs_tab_widget = QWidget()
        logs_tab_layout = QVBoxLayout()
        self.logs_label = QLabel("")
        self.logs_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        scrollable_view = QScrollArea()
        scrollable_view.setWidgetResizable(True)
        scrollable_view.setWidget(self.logs_label)
        logs_tab_layout.addWidget(scrollable_view)
        logs_tab_widget.setLayout(logs_tab_layout)
        logs_tab_widget.setStyleSheet("background-color: #FFFFFF")
        return logs_tab_widget

    def settings_tab(self):
        settings_tab_widget = QWidget()
        settings_tab_layout = QVBoxLayout()
        self.settings_label = QLabel("")
        self.settings_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        scrollable_view = QScrollArea()
        scrollable_view.setWidgetResizable(True)
        scrollable_view.setWidget(self.settings_label)
        settings_tab_layout.addWidget(scrollable_view)
        settings_tab_widget.setLayout(settings_tab_layout)
        settings_tab_widget.setStyleSheet("background-color: #FFFFFF")
        self.settings_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        return settings_tab_widget
    
    def transcript_tab(self, transcript_name, blast_results):
        transcript_tab_widget = QWidget()
        transcript_tab_layout = QVBoxLayout()
        inner_transcript_layout = QVBoxLayout()
        scrollable_view = QScrollArea()
        scrollable_view.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scrollable_view.setWidgetResizable(True)
        transcript_tab_layout.addWidget(scrollable_view)
        transcript_inner_widget = QWidget()
        inner_transcript_layout.addWidget(QLabel("------------------------------------------- \n"
                                                 "Transcript Information \n"
                                                 "------------------------------------------- \n"))
        inner_transcript_layout.addWidget(TranscriptResultView(blast_results[0]))
        inner_transcript_layout.addWidget(QLabel("------------------------------------------- \n"
                                                 "BLAST RESULTS \n"
                                                 "------------------------------------------- \n"))
        for result in blast_results:
            inner_transcript_layout.addWidget(BlastResultView(result))

        transcript_inner_widget.setLayout(inner_transcript_layout)
        scrollable_view.setWidget(transcript_inner_widget)
        transcript_tab_widget.setLayout(transcript_tab_layout)
        transcript_tab_widget.setStyleSheet("background-color: #FFFFFF")

        self.tabs_widget.addTab(transcript_tab_widget, transcript_name)