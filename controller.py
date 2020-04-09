import argparse

from helpers.accession_handler import *
from helpers.filter_handler import *
from helpers.setup_helper import handle_setup

from helpers.blast_helper import *

# define arguments for the script.
parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('--name', metavar='name', type=str, nargs='+',
                    help='(required) a name for your test run')
parser.add_argument('--quiet', metavar='quiet', type=str, nargs='+',
                    help='y to only log to log file, n (default) to log to log file and print to console.', default="n")
parser.add_argument('--test_data', metavar='test_data', type=str, nargs='+',
                    help='y (default) to use sample data, n to download new files.', default="y")


def handle_transcripts_post_filtering(global_data):
    if len(global_data.valid_transcripts) == 0:
        global_data.logger.log("No valid transcripts are left after filtering")
        exit()
    else:
        global_data.logger.log(f"Valid transcripts after filtering {global_data.valid_transcripts}")


global_data, global_settings = handle_setup(parser=parser)
# run all of the functions from the helper scripts
global_data = get_gene(global_data)
global_data = separate_gene_transcripts(global_data)
global_data = get_promoter_sequence(global_data)
global_data = filter_transcripts(global_data)
handle_transcripts_post_filtering(global_data)

blast_local = False
if blast_local:
    global_data = locally_handle_blasting_transcripts(global_data)
    global_data = get_blasted_transcript_information(global_data)
else:
    global_data = remotely_handle_blasting_transcripts(global_data)

export_blasted_results(global_data)
global_data.export_as_json()
global_settings.export_as_json()
global_data.logger.log("Completed running all scripts")
