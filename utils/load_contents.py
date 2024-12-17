
import os

def load_contents(path_to_file):
    if not exists(path_to_file):
        raise FileNotFoundError(f"The file at '{path_to_file}' does not exist.")
    try:
        with open(path_to_file, 'r', encoding='utf-8') as file:
            content = file.read()
        return content
    except Exception as e:
        raise IOError(f"An error occurred while reading the file: {e}")

