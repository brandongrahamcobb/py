from typing import Any, Dict

global logger

import load_yaml
import logging
import logging.handlers


def setup_logging(config: Dict[str, Any]) -> None:
    config = load_yaml(helpers.hapth_config_yaml)
    logging_level = config['logging_level'].upper()
    logging.basicConfig(level=getattr(logging, logging_level))
    if not exists(dirname(path_log)):
        makedirs(dirname(path_log))
    file_handler = logging.FileHandler(path_log)
    file_handler.setLevel(logging_level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger = logging.getLogger()
    logger.setLevel(logging_level)
    logger.addHandler(file_handler)
