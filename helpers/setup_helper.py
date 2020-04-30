import subprocess

from helpers.log_handler import Logger
from helpers.input_data_handler import *
import gzip
import shutil
import urllib.request as request
from contextlib import closing


def create_test_folder(folder_name):
    # make sure the folder name is valid
    if folder_name is None or folder_name == "":
        print("Invalid folder name")
        return None

    folder_path = os.path.join(os.getcwd(), f"{FOLDER_PREFIX}{folder_name}")
    try:
        os.mkdir(folder_path)
    except Exception as e:
        print("FOLDER ALREADY EXISTS")
    return folder_path


def _handle_run_folder_creation(args):
    final_folder_path = create_test_folder(str(args.name[0]).replace(' ', '_'))
    # exit if folder is not made.
    if not final_folder_path:
        exit()
    return final_folder_path


def _export_default_global_settings(folder_path):
    if os.path.exists(os.path.join(folder_path, "run_settings.json")):
        global_settings = GlobalApplicationSettings(folder_path)
    else:
        global_settings = GlobalApplicationSettings(folder_path)
        global_settings.add_to_settings(FILTER_TRANSCRIPT_INCOMPLETE, True)
        global_settings.add_to_settings(FILTER_TRANSCRIPT_PROTEIN_CODING, True)
        global_settings.add_to_settings(FILTER_TRANSCRIPT_START_CODON, True)
        global_settings.add_to_settings(FILTER_TRANSCRIPT_STOP_CODON, True)
        global_settings.export_as_json()
    return global_settings


def handle_setup(gui=True, folder_path=None, parser=None):
    if gui:
        final_folder_path = create_test_folder(folder_path)
        logger = Logger(final_folder_path, True)
        global_data = GlobalApplicationData(final_folder_path, logger, "LOCAL", is_gui=gui, is_writing=True)
        global_settings = GlobalApplicationSettings(final_folder_path)
        return global_data, global_settings, final_folder_path
    else:
        args = parser.parse_args()
        final_folder_path = _handle_run_folder_creation(args)
        logger = Logger(final_folder_path, args.quiet[0])
        input_global_data = GlobalApplicationData(final_folder_path, logger, args.blast_type, is_writing=True,
                                                  args=args)
        global_settings = _export_default_global_settings(final_folder_path)
        return input_global_data, global_settings


def handle_download_cdna_locally():
    if not os.path.exists(os.path.join(os.getcwd(), "zebrafish_cdna_from_genomic.fna.gz")):
        with closing(request.urlopen('https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Danio_rerio/all_assembly_versions/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_cds_from_genomic.fna.gz')) as r:
            with open(os.path.join(os.getcwd(), "zebrafish_cdna_from_genomic.fna.gz"), 'wb') as f:
                shutil.copyfileobj(r, f)

    if not os.path.exists(os.path.join(os.getcwd(), "zebrafish_cdna_from_genomic.fasta")):
        with gzip.open(os.path.join(os.getcwd(), "zebrafish_cdna_from_genomic.fna.gz"), "rb") as f_in:
            with open(os.path.join(os.getcwd(), "zebrafish_cdna_from_genomic.fasta"), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    os.system(f"makeblastdb -in zebrafish_cdna_from_genomic.fasta -parse_seqids -title zebrafish_blastdb"
              f" -dbtype nucl")

