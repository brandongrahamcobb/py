
#    Lucy, a discord.py bot, is an open-source package containing four cogs and a helper file.
#    Copyright (C) 2024  Cobb, Brandon Graham

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

from os import getenv, makedirs
from os.path import abspath, dirname, exists, expanduser, join

import asyncio
import datetime as dt
import os
import shutil
import subprocess
import sys
import tarfile

class Project:

    def __init__(self):

        self.COBBBRANDONGRAHAM = expanduser('~')

        self.CURRENT_DATE = dt.datetime.now().strftime("%d%m%y")
        self.PATH_DOCUMENTS_MD_README = join(self.COBBBRANDONGRAHAM, 'Documents', 'md', f'CobbBrandonGraham_README_{self.CURRENT_DATE}.md')
        self.PATH_DOCUMENTS_PY = abspath(__file__)
        self.PATH_DOCUMENTS_PY_ADMIN_COG = join(self.COBBBRANDONGRAHAM, 'Documents', 'py', 'bot', 'cogs', f'CobbBrandonGraham_admin_cog_{self.CURRENT_DATE}.py')
        self.PATH_DOCUMENTS_PY_AI_COG = join(self.COBBBRANDONGRAHAM, 'Documents', 'py', 'bot', 'cogs', f'CobbBrandonGraham_ai_cog_{self.CURRENT_DATE}.py')
        self.PATH_DOCUMENTS_PY_COGS_INIT = join(self.COBBBRANDONGRAHAM, 'Documents', 'py', 'bot', 'cogs', f'CobbBrandonGraham__init__{self.CURRENT_DATE}.py')
        self.PATH_DOCUMENTS_PY_INIT = join(self.COBBBRANDONGRAHAM, 'Documents', 'py', 'bot', f'CobbBrandonGraham__init__{self.CURRENT_DATE}.py')
        self.PATH_DOCUMENTS_PY_HELPERS = join(self.COBBBRANDONGRAHAM, 'Documents', 'py', 'bot', 'utils', f'CobbBrandonGraham_helpers_{self.CURRENT_DATE}.py')
        self.PATH_DOCUMENTS_PY_MAIN = join(self.COBBBRANDONGRAHAM, 'Documents', 'py', 'bot', f'CobbBrandonGraham_main_{self.CURRENT_DATE}.py')
        self.PATH_DOCUMENTS_PY_MY_COG = join(self.COBBBRANDONGRAHAM, 'Documents', 'py', 'bot', 'cogs', f'CobbBrandonGraham_my_cog_{self.CURRENT_DATE}.py')
        self.PATH_DOCUMENTS_PY_SETUP = join(self.COBBBRANDONGRAHAM, 'Documents', 'py', f'CobbBrandonGraham_setup_{self.CURRENT_DATE}.py')
        self.PATH_DOCUMENTS_PY_USER_COG = join(self.COBBBRANDONGRAHAM, 'Documents', 'py', 'bot', 'cogs', f'CobbBrandonGraham_user_cog_{self.CURRENT_DATE}.py')
        self.PATH_DOCUMENTS_PY_UTILS_INIT = join(self.COBBBRANDONGRAHAM, 'Documents', 'py', 'bot', 'utils', f'CobbBrandonGraham__init__{self.CURRENT_DATE}.py')
        self.PATH_DOCUMENTS_TXT_LICENSE = join(self.COBBBRANDONGRAHAM, 'Documents', 'txt', f'GNU_LICENSE_{self.CURRENT_DATE}.txt')
        self.PATH_DOCUMENTS_TXT_REQUIREMENTS = join(self.COBBBRANDONGRAHAM, 'Documents', 'txt', f'CobbBrandonGraham_requirements_{self.CURRENT_DATE}.txt')

        self.PATH_SRC_PKGBUILD = join(self.COBBBRANDONGRAHAM, 'Documents', 'PKGBUILD')
        self.PATH_SRC_LICENSE = join(self.COBBBRANDONGRAHAM, 'Documents', 'src', f'lucy', 'LICENSE')
        self.PATH_SRC_PY_ADMIN_COG = join(self.COBBBRANDONGRAHAM, 'Documents', 'src', f'lucy', 'bot', 'cogs', 'admin_cog.py')
        self.PATH_SRC_PY_AI_COG = join(self.COBBBRANDONGRAHAM, 'Documents', 'src', f'lucy', 'bot', 'cogs', 'ai_cog.py')
        self.PATH_SRC_PY_COGS_INIT = join(self.COBBBRANDONGRAHAM, 'Documents', 'src', f'lucy', 'bot', 'cogs', '__init__.py')
        self.PATH_SRC_PY_HELPERS = join(self.COBBBRANDONGRAHAM, 'Documents', 'src', f'lucy', 'bot', 'utils', 'helpers.py')
        self.PATH_SRC_PY_INIT = join(self.COBBBRANDONGRAHAM, 'Documents', 'src', f'lucy', 'bot', '__init__.py')
        self.PATH_SRC_PY_MAIN = join(self.COBBBRANDONGRAHAM, 'Documents', 'src', f'lucy', 'bot', 'main.py')
        self.PATH_SRC_PY_MY_COG = join(self.COBBBRANDONGRAHAM, 'Documents', 'src', f'lucy', 'bot', 'cogs', 'my_cog.py')
        self.PATH_SRC_PY_SETUP = join(self.COBBBRANDONGRAHAM, 'Documents', 'src', f'lucy', 'setup.py')
        self.PATH_SRC_PY_USER_COG = join(self.COBBBRANDONGRAHAM, 'Documents', 'src', f'lucy', 'bot', 'cogs', 'user_cog.py')
        self.PATH_SRC_PY_UTILS_INIT = join(self.COBBBRANDONGRAHAM, 'Documents', 'src', f'lucy', 'bot', 'utils', '__init__.py')
        self.PATH_SRC_TXT_REQUIREMENTS = join(self.COBBBRANDONGRAHAM, 'Documents', 'src', f'lucy', 'requirements.txt')
        self.PATH_SRC_README = join(self.COBBBRANDONGRAHAM, 'Documents', 'src', f'lucy', 'README.md')

def run_makepkg():
    try:
        result = subprocess.run(['makepkg', '-si'], check=True, text=True, capture_output=True)
        print('Command Output:', result.stdout)
        print('Command Error:', result.stderr)
    except subprocess.CalledProcessError as e:
        print('An error occurred while running makepkg:')
        print('Error Code:', e.returncode)
        print('Error Output:', e.stderr)

def verify_copy(input, src):

    if exists(src):

        print(f'Copied {input} to {src}')

    else:

        print(f'Failed to copy {input}')

if __name__ == '__main__':

    self = Project()
    makedirs(dirname(self.PATH_SRC_README), exist_ok=True)
    makedirs(dirname(self.PATH_SRC_PY_ADMIN_COG), exist_ok=True)
    makedirs(dirname(self.PATH_SRC_PY_AI_COG), exist_ok=True)
    makedirs(dirname(self.PATH_SRC_PY_COGS_INIT), exist_ok=True)
    makedirs(dirname(self.PATH_SRC_PY_HELPERS),	 exist_ok=True)
    makedirs(dirname(self.PATH_SRC_PY_INIT), exist_ok=True)
    makedirs(dirname(self.PATH_SRC_PY_MAIN), exist_ok=True)
    makedirs(dirname(self.PATH_SRC_PY_MY_COG), exist_ok=True)
    makedirs(dirname(self.PATH_SRC_PY_SETUP), exist_ok=True)
    makedirs(dirname(self.PATH_SRC_PY_USER_COG), exist_ok=True)
    makedirs(dirname(self.PATH_SRC_PY_UTILS_INIT), exist_ok=True)
    makedirs(dirname(self.PATH_SRC_LICENSE), exist_ok=True)
    makedirs(dirname(self.PATH_SRC_TXT_REQUIREMENTS), exist_ok=True)

    shutil.copy2(self.PATH_DOCUMENTS_MD_README, self.PATH_SRC_README)
    shutil.copy2(self.PATH_DOCUMENTS_PY_ADMIN_COG, self.PATH_SRC_PY_ADMIN_COG)
    shutil.copy2(self.PATH_DOCUMENTS_PY_AI_COG, self.PATH_SRC_PY_AI_COG)
    shutil.copy2(self.PATH_DOCUMENTS_PY_COGS_INIT, self.PATH_SRC_PY_COGS_INIT)
    shutil.copy2(self.PATH_DOCUMENTS_PY_HELPERS, self.PATH_SRC_PY_HELPERS)
    shutil.copy2(self.PATH_DOCUMENTS_PY_INIT, self.PATH_SRC_PY_INIT)
    shutil.copy2(self.PATH_DOCUMENTS_PY_MAIN, self.PATH_SRC_PY_MAIN)
    shutil.copy2(self.PATH_DOCUMENTS_PY_MY_COG, self.PATH_SRC_PY_MY_COG)
    shutil.copy2(self.PATH_DOCUMENTS_PY_SETUP, self.PATH_SRC_PY_SETUP)
    shutil.copy2(self.PATH_DOCUMENTS_PY_USER_COG, self.PATH_SRC_PY_USER_COG)
    shutil.copy2(self.PATH_DOCUMENTS_PY_UTILS_INIT, self.PATH_SRC_PY_UTILS_INIT)
    shutil.copy2(self.PATH_DOCUMENTS_TXT_LICENSE, self.PATH_SRC_LICENSE)
    shutil.copy2(self.PATH_DOCUMENTS_TXT_REQUIREMENTS, self.PATH_SRC_TXT_REQUIREMENTS)

    verify_copy(self.PATH_DOCUMENTS_MD_README, self.PATH_SRC_README)
    verify_copy(self.PATH_DOCUMENTS_PY_ADMIN_COG, self.PATH_SRC_PY_ADMIN_COG)
    verify_copy(self.PATH_DOCUMENTS_PY_AI_COG, self.PATH_SRC_PY_AI_COG)
    verify_copy(self.PATH_DOCUMENTS_PY_COGS_INIT, self.PATH_SRC_PY_COGS_INIT)
    verify_copy(self.PATH_DOCUMENTS_PY_HELPERS, self.PATH_SRC_PY_HELPERS)
    verify_copy(self.PATH_DOCUMENTS_PY_INIT, self.PATH_SRC_PY_INIT)
    verify_copy(self.PATH_DOCUMENTS_PY_MAIN, self.PATH_SRC_PY_MAIN)
    verify_copy(self.PATH_DOCUMENTS_PY_MY_COG, self.PATH_SRC_PY_MY_COG)
    verify_copy(self.PATH_DOCUMENTS_PY_SETUP, self.PATH_SRC_PY_SETUP)
    verify_copy(self.PATH_DOCUMENTS_PY_USER_COG, self.PATH_SRC_PY_USER_COG)
    verify_copy(self.PATH_DOCUMENTS_PY_UTILS_INIT, self.PATH_SRC_PY_UTILS_INIT)
    verify_copy(self.PATH_DOCUMENTS_TXT_LICENSE, self.PATH_SRC_LICENSE)
    verify_copy(self.PATH_DOCUMENTS_TXT_REQUIREMENTS, self.PATH_SRC_TXT_REQUIREMENTS)

    with tarfile.open(f'CobbBrandonGraham_lucy_{self.CURRENT_DATE}.tar.gz', 'w:gz') as tar:

        tar.add(self.PATH_SRC_README)
        tar.add(self.PATH_SRC_PKGBUILD)
        tar.add(self.PATH_SRC_PY_ADMIN_COG)
        tar.add(self.PATH_SRC_PY_AI_COG)
        tar.add(self.PATH_SRC_PY_COGS_INIT)
        tar.add(self.PATH_SRC_PY_HELPERS)
        tar.add(self.PATH_SRC_PY_INIT)
        tar.add(self.PATH_SRC_PY_MAIN)
        tar.add(self.PATH_SRC_PY_MY_COG)
        tar.add(self.PATH_SRC_PY_SETUP)
        tar.add(self.PATH_SRC_PY_USER_COG)
        tar.add(self.PATH_SRC_PY_UTILS_INIT)
        tar.add(self.PATH_SRC_LICENSE)
        tar.add(self.PATH_SRC_TXT_REQUIREMENTS)
        tar.close()

#    run_makepkg()    #package_path = dirname(abspath(__file__))
   # if package_path not in sys.path:
  #      sys.path.insert(0, package_path)

 #   import bot.main
#    asyncio.run(bot.main.main())
