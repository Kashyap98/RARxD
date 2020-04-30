import os
import json
import jsonpickle

from helpers.constants import *


class GlobalApplicationData(object):

    def __init__(self, folder_path, logger, blast_type, is_gui=False, is_writing=True, args=None):
        self.folder_path = folder_path
        self.handle_gui_command_line_input(args, is_gui)
        self.input_keys = self.input_data.keys()
        self.input_data = self.input_data
        self.logger = logger
        self.target_sequence = self.handle_input_data(TARGET_SEQUENCE)
        self.input_promoter = self.handle_input_data(PROMOTER)
        self.promoter_start = self.handle_input_data(PROMOTER_START)
        self.promoter_end = self.handle_input_data(PROMOTER_END)
        self.ensembl_gene = None
        self.transcripts = {}
        self.valid_transcripts = []
        self.invalid_transcripts = []
        self.promoter = None
        self.transcript_results = []
        if is_gui:
            self.blast_type = blast_type
        else:
            self.blast_type = self.handle_input_data(BLAST_TYPE)
        self.blast_local = self.handle_blast_type(self.blast_type)
        self.logger.log(f"Input data found: {self.input_data}")

    def handle_gui_command_line_input(self, args, is_gui):
        if is_gui:
            try:
                self.input_data = self.load_json(self.folder_path)
            except Exception as e:
                self.input_data = {}
        else:
            if args.gene_id[0] is None:
                self.logger.log("Missing Gene ID. Required to run this program")
                exit()
            if args.promoter_id[0] is None:
                self.logger.log("Missing Promoter ID. Required to run this program")
                exit()

            self.input_data = {TARGET_SEQUENCE: args.gene_id[0], PROMOTER: args.promoter_id[0],
                               BLAST_TYPE: args.blast_type[0]}

    def handle_blast_type(self, blast_type):
        if blast_type == "L":
            return True
        else:
            return False

    def handle_input_data(self, input_key):
        requested_data = None
        if input_key in self.input_keys:
            requested_data = self.input_data[input_key]
        return requested_data

    def invalidate_transcript(self, transcript_name):
        self.valid_transcripts.remove(transcript_name)
        self.invalid_transcripts.append(transcript_name)

    def load_json(self, folder_path):
        with open(os.path.join(folder_path, f"{INPUT_DATA_FILE_NAME}"), "r") as sample_file:
            return json.loads(sample_file.read())

    def export_as_json(self):
        with open(os.path.join(self.folder_path, f"{INPUT_DATA_FILE_NAME}"), "w") as output_file:
            output_file.write(jsonpickle.encode(self, keys=True))


class GlobalApplicationSettings(object):

    def __init__(self, folder_path):
        self.folder_path = folder_path
        self.file_path = os.path.join(self.folder_path, f"{SETTINGS_FILE_NAME}")
        self.input_settings = {}
        self.load_json()

    def add_to_settings(self, setting_name, setting_value):
        self.input_settings[setting_name] = setting_value

    def load_json(self):
        # Replace path when ready (self.global_data.folder_path, f"{SETTINGS_FILE_NAME}")
        if os.path.exists(self.file_path):
            with open(self.file_path, "r") as input_file:
                file_data = json.loads(input_file.read())
                if file_data != {}:
                    self.input_settings = file_data
                # self.folder_path = self.input_settings["folder_path"]

    def export_as_json(self):
        with open(self.file_path, "w") as output_file:
            output_file.write(jsonpickle.encode(self.input_settings, keys=True))
