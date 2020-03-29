import os
import json

from helpers.constants import *


class GlobalApplicationData:

    def __init__(self, folder_path, logger, input_data):
        self.input_keys = input_data.keys()
        self.input_data = input_data
        self.folder_path = folder_path
        self.logger = logger
        self.target_sequence = self.handle_input_data(TARGET_SEQUENCE)
        self.promoter = self.handle_input_data(PROMOTER)
        self.promoter_start = self.handle_input_data(PROMOTER_START)
        self.promoter_end = self.handle_input_data(PROMOTER_END)
        self.ensembl_gene = None
        self.transcripts = {}
        self.valid_transcripts = []
        self.invalid_transcripts = []

        self.logger.log(f"Input data found: {input_data}")
        
    def handle_input_data(self, input_key):
        requested_data = None
        if input_key in self.input_keys:
            requested_data = self.input_data[input_key]
        return requested_data

    def invalidate_transcript(self, transcript_name):
        self.valid_transcripts.remove(transcript_name)
        self.invalid_transcripts.append(transcript_name)


def handle_json_data_input(folder_path, logger):
    with open(os.path.join(os.getcwd(), "sample_data.txt"), "r") as sample_file:
        input_data = json.loads(sample_file.read())
        global_data = GlobalApplicationData(folder_path, logger, input_data)
        return global_data
