import os


def get_quiet_mode(quiet_arg):
    if quiet_arg == "y":
        return True
    else:
        return False


class Logger:
    # init logger object
    def __init__(self, folder_path, quiet_arg):
        self.folder_path = folder_path
        self.log_file_name = os.path.basename(folder_path[12:])
        self.file_path = os.path.join(self.folder_path, f"{self.log_file_name}_log.txt")
        self.quiet_mode = get_quiet_mode(quiet_arg)

        self.create_log_file()

    # responsible for all logging, check if outputting to log and terminal
    def log(self, input_log):
        with open(self.file_path, "a") as log_file:
            log_file.write(f"{input_log}\n")
            if not self.quiet_mode:
                print(str(input_log))

    # initial creation of log file
    def create_log_file(self):
        self.log("Creating Log File")

    # set quiet_mode arg
