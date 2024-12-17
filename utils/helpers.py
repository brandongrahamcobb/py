''' helpers.py  The purpose of this program is to provide generic directiories.
    Copyright (C) 2024 github.com/brandongrahamcobb

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''
from os.path import dirname, abspath, expanduser, join

# Base and Home Paths
DIR_BASE = dirname(abspath(__file__))
PATH_HOME = expanduser('~')

# Config and Preferences
PATH_CONFIG_YAML = join(PATH_HOME, '.config', 'spawd', 'config.yaml')
PATH_CONFIG_JSON = join(PATH_HOME, '.config', 'spawd', 'config.json')

# Logs
PATH_LOG = join(PATH_HOME, '.log', 'spawd', 'discord.log')

# Script Paths
PATH_ADJUST_HUE_AND_SATURATION = 'adjust_hue_and_saturation.py'
PATH_ARPP = 'api_request_parallel_processor.py'
PATH_BOT = 'bot.py'
PATH_COMBINE = 'combine.py'
PATH_CREATE_BATCH_COMPLETION = 'create_batch_completion.py'
PATH_CREATE_COMPLETION = 'create_completion.py'
PATH_CREATE_HTTPS_COMPLETION = 'create_https_completion.py'
PATH_CREATE_MODERATION = 'create_moderation.py'
PATH_DRAW_FINGERPRINT = 'draw_fingerprint.py'
PATH_DRAW_WATERMARKED_MOLECULE = 'draw_watermarked_molecule.py'
PATH_FORMAT_ERROR_CHECK = 'format_error_check.py'
PATH_GET_MOLECULE_NAME = 'get_molecule_name.py'
PATH_GET_MOL = 'get_mol.py'
PATH_GET_PROXIMITY = 'get_proximity.py'
PATH_GET_SCRIPTURE = 'get_scripture.py'
PATH_GOOGLE = 'google.py'
PATH_GSRS = 'gsrs.py'
PATH_LOAD_CONTENTS = 'load_contents.py'
PATH_LOAD_YAML = 'load_yaml.py'
PATH_SETUP_LOGGING = 'setup_logging.py'
PATH_STABLE_CASCADE = 'stable_cascade.py'
PATH_UNIQUE_PAIRS = 'unique_pairs.py'
PATH_VYRTUOUS = 'vyrtuous.py'
