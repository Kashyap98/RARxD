from helpers.log_handler import Logger
from helpers.input_data_handler import *
from helpers.blast_helper import *


def create_test_folder(folder_name, gui=False):
    # make sure the folder name is valid
    if folder_name is None or folder_name == "":
        print("Invalid folder name")
        return None

    if gui:
        folder_path = os.path.join(os.getcwd(), "..", f"{FOLDER_PREFIX}{folder_name}")
    else:
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


def handle_setup(gui=False, folder_path=None, parser=None):
    if gui:
        final_folder_path = create_test_folder(folder_path, gui)
        logger = Logger(final_folder_path, True)
        return final_folder_path
    else:
        args = parser.parse_args()
        final_folder_path = create_test_folder(str(args.name[0]).replace(' ', '_'))
        # exit if folder is not made.
        if not final_folder_path:
            exit()
        logger = Logger(final_folder_path, args.quiet[0])
        input_global_data = GlobalApplicationData(final_folder_path, logger, is_writing=True)
        global_settings = GlobalApplicationSettings(final_folder_path, is_writing=True)
        return input_global_data, global_settings
