import os
import shutil
import sys
import glob
import argparse

from helpers.constants import *

YES = ["yes", "y"]
NO = ["no", "n"]


def get_user_choice():
    while True:
        choice = input().lower()
        if choice in YES:
            return True
        elif choice in NO:
            return False
        else:
            sys.stdout.write("Please respond with 'yes' or 'no'")


print("Are you sure you want to delete all miniProject test directories? This cannot be undone.")
print("[y/n]")


def remove_folder(folder_path):
    folder_deleted = False
    if os.path.isdir(folder_path):
        shutil.rmtree(folder_path)
        if not os.path.exists(folder_path):
            folder_deleted = True
    return folder_deleted


def delete_folders(all_folders):
    print("Deleting folders.")
    delete_count = 0
    for folder in all_folders:
        if remove_folder(folder):
            delete_count += 1
    print(f"Deleted {delete_count} folders.")


def run_cleaner(ask_user_choice=True):
    all_folders = glob.glob(os.path.join(os.getcwd(), f"{FOLDER_PREFIX}*"))
    if ask_user_choice:
        if get_user_choice():
            delete_folders(all_folders)
        else:
            print("Not deleting any folders.")
    else:
        delete_folders(all_folders)


def main():
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('--force', '-force', '-f', action='store_true',
                        help='(optional) force folder delete, do not ask user for choice')
    args = parser.parse_args()
    if args.force is True:
        run_cleaner(ask_user_choice=False)
    else:
        run_cleaner()


if __name__ == "__main__":
    main()
