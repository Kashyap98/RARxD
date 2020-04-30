import argparse
import sys

from PySide2.QtWidgets import QApplication

from gui.progress_bar_util import ProgressBarThread
from helpers.accession_handler import *
from helpers.filter_handler import *
from helpers.setup_helper import handle_setup
from helpers.output_handler import export_genbank_files
from helpers.blast_helper import *

# define arguments for the script.
parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('--name', metavar='name', type=str, nargs='+',
                    help='(required) a name for your test run')
parser.add_argument('--quiet', metavar='quiet', type=str, nargs='+',
                    help='y to only log to log file, n (default) to log to log file and print to console.', default="n")
parser.add_argument('--blast_type', metavar='blast_type', type=str, nargs='+',
                    help='LOCAL (default) to local blastn database, REMOTE to remotely BLAST NCBI server.',
                    default="LOCAL")
parser.add_argument('--gui_mode', metavar='gui_mode', type=str, nargs='+',
                    help='n (default) when not run from gui, y when running from gui.',
                    default="n")
parser.add_argument('--gene_id', metavar='gene_id', type=str, nargs='+',
                    help='Ensembl Accession ENSDARGXXXXXXXXXXX:')
parser.add_argument('--promoter_id', metavar='promoter_id', type=str, nargs='+',
                    help='Promoter accession from NCBI to prepend to the transcript.')


def handle_transcripts_post_filtering(global_data):
    if len(global_data.valid_transcripts) == 0:
        global_data.logger.log("No valid transcripts are left after filtering")
        exit()
    else:
        global_data.logger.log(f"Valid transcripts after filtering {global_data.valid_transcripts}")


app = QApplication(sys.argv)
progress_bar = ProgressBarThread(app=app, max_completed=10)

global_data, global_settings = handle_setup(parser=parser, gui=False)
progress_bar.increment_progress("Finished set up tasks")
# run all of the functions from the helper scripts
global_data = get_gene(global_data)
progress_bar.increment_progress("Retrieved gene from Ensembl")
global_data = separate_gene_transcripts(global_data)
progress_bar.increment_progress("Separated gene transcripts")
global_data = get_promoter_sequence(global_data)
progress_bar.increment_progress("Retrieved promoter from NCBI")
global_data = filter_transcripts(global_data, global_settings)
progress_bar.increment_progress("Filtered transcripts")
handle_transcripts_post_filtering(global_data)
progress_bar.increment_progress("Post filter transcript validation")

if global_data.blast_local:
    global_data = locally_handle_blasting_transcripts(global_data)
    global_data = get_blasted_transcript_information(global_data)
else:
    global_data = remotely_handle_blasting_transcripts(global_data)
progress_bar.increment_progress("Completed BLAST analysis")
export_blasted_results(global_data)
progress_bar.increment_progress("Exported BLAST results")
export_genbank_files(global_data)
progress_bar.increment_progress("Exported Genbank results")
global_data.logger.log("Completed running all scripts")
progress_bar.increment_progress("Done")

