
import os
import sys

def clear_screen():
    if sys.platform == "win32":
        os.system('cls')
    else:
        os.system('clear')
