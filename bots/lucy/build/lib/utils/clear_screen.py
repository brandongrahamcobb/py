''' clear_screen.py  The purpose of this program is clear screen on terminals either in Windows or Linux functionality from cd ../.
    Copyright (C) 2024  github.com/brandongrahamcobb

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

from utils.setup_logging import logger

import os
import sys

def clear_screen():
    logger.info('Clearing screen.')
    if sys.platform == 'win32':
        os.system('cls')
    else:
        os.system('clear')