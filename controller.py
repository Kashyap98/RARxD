import argparse
import os

from helpers.log_handler import Logger

# define arguments for the script.
parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('--name', metavar='name', type=str, nargs='+',
                    help='(required) a name for your test run')
parser.add_argument('--quiet', metavar='quiet', type=str, nargs='+',
                    help='y to only log to log file, n (default) to log to log file and print to console.', default="n")
parser.add_argument('--test_data', metavar='test_data', type=str, nargs='+',
                    help='y (default) to use sample data, n to download new files.', default="y")


def create_test_folder(args):
    folder_name = str(args.name[0]).replace(' ', '_')
    # make sure the folder name is valid
    if folder_name is None or folder_name == "":
        print("Invalid folder name")
        return None
    folder_path = os.path.join(os.getcwd(), f"miniProject_{folder_name}")

    # check if folder exists, if not make it.
    if os.path.exists(folder_path):
        print("Folder name already exits. Please run cleaner.py")
    else:
        os.mkdir(folder_path)
        return folder_path


def arg_get_files(input_arg):
    # Determine if downloading files using wget
    if input_arg == "y":
        return True
    else:
        return False


args = parser.parse_args()
folder_path = create_test_folder(args)
# exit if folder is not made.
if not folder_path:
    exit()

# run the helpers according to the instructions.
logger = Logger(folder_path, args.quiet[0])
logger.log("Completed running all helpers")