from typing import Any, Dict

import yaml

def increment_version(config: Dict[str, Any], path_config_yaml):
    current_version = config['version']
    major, minor, patch = map(int, current_version.split('.'))
    patch += 1
    if patch >= 10:
        patch = 0
        minor += 1
    if minor >= 10:
        minor = 0
        major += 1
    new_version = f'{major}.{minor}.{patch}'
    config['version'] = new_version
    with open(path_config_yaml, 'w') as file:
        yaml.dump(config, file)
