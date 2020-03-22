import os
import shutil
import sys
import glob

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

if get_user_choice():
    print("Deleting folders.")
    delete_count = 0
    for folder in glob.glob(os.path.join(os.getcwd(), "run_*")):
        if os.path.isdir(folder):
            shutil.rmtree(folder)
            delete_count += 1
    print(f"Deleted {delete_count} folders.")

else:
    print("Not deleting any folders.")
