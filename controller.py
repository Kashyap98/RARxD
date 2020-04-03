import argparse

from helpers.log_handler import Logger
from helpers.accession_handler import *
from helpers.input_data_handler import *
from helpers.filter_handler import *
from helpers.blast_helper import *

# define arguments for the script.
parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('--name', metavar='name', type=str, nargs='+',
                    help='(required) a name for your test run')
parser.add_argument('--quiet', metavar='quiet', type=str, nargs='+',
                    help='y to only log to log file, n (default) to log to log file and print to console.', default="n")
parser.add_argument('--test_data', metavar='test_data', type=str, nargs='+',
                    help='y (default) to use sample data, n to download new files.', default="y")


def create_test_folder(input_args):
    folder_name = str(input_args.name[0]).replace(' ', '_')
    # make sure the folder name is valid
    if folder_name is None or folder_name == "":
        print("Invalid folder name")
        return None
    folder_path = os.path.join(os.getcwd(), f"{FOLDER_PREFIX}{folder_name}")

    # # check if folder exists, if not make it.
    if os.path.exists(folder_path):
        print("Folder name already exits. Please run cleaner.py")
    else:
        # TODO REMOVE LATER
        try:
            os.mkdir(folder_path)
        except Exception as e:
            print("FOLDER ALREADY EXISTS")
    return folder_path


def handle_setup():
    args = parser.parse_args()
    final_folder_path = create_test_folder(args)
    # exit if folder is not made.
    if not final_folder_path:
        exit()
    logger = Logger(final_folder_path, args.quiet[0])
    input_global_data = GlobalApplicationData(final_folder_path, logger, is_writing=False)
    global_settings = GlobalApplicationSettings(final_folder_path, is_writing=False)
    return input_global_data, global_settings


def handle_transcripts_post_filtering(global_data):
    if len(global_data.valid_transcripts) == 0:
        global_data.logger.log("No valid transcripts are left after filtering")
        exit()
    else:
        global_data.logger.log(f"Valid transcripts after filtering {global_data.valid_transcripts}")


global_data, global_settings = handle_setup()
# run all of the functions from the helper scripts
global_data = get_gene(global_data)
global_data = separate_gene_transcripts(global_data)
global_data = get_promoter_sequence(global_data)
global_data = filter_transcripts(global_data)
handle_transcripts_post_filtering(global_data)

blast_local = True
if blast_local:
    global_data = locally_handle_blasting_transcripts(global_data)
    global_data = get_blasted_transcript_information(global_data)
else:
    global_data = remotely_handle_blasting_transcripts(global_data)

global_data.logger.log("Completed running all scripts")
