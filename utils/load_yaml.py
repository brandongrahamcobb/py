import os
import traceback
import yaml

def load_yaml(path_to_file):
    try:
        if not os.path.exists(path_to_file):
            return {}
        with open(path_to_file, 'r') as f:
            return yaml.safe_load(f) or {}
    except Exception as e:
        print(traceback.format_exc())
