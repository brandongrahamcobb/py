'''
    Copyright (C) 2024 https://github.com/brandongrahamcobb

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <https://www.gnu.org/licenses/>.

'''

# ./run.py

import json
import os
import requests
import subprocess
import sys

requirements = ["discord.py", "pubchempy", "rdkit", "pillow", "requests"]

dir = os.path.dirname(os.path.abspath(__file__))
updir = os.path.join(dir, '..')
dir_json = os.path.join(updir, 'json')
dir_log = os.path.join(updir, 'log')
dir_txt = os.path.join(updir, 'txt')
dir_venv = os.path.join(updir, 'activate')
file_json = os.path.join(dir_json, 'config.json')
file_py = os.path.join(dir, 'advanced_startup.py')
system = 'Linux' if os.name != 'nt' else 'Windows'
py = os.path.join(dir_venv, 'bin', 'python') if system == 'Linux' else sys.executable

if __name__ == '__main__':
    print(f'System: {system}')
    print(f'Location: {dir}')
    print(f'Python Location: {sys.executable}')
    print(f'Python Version: {sys.version}')
    print('Updating packages whether you like it or not...')
    subprocess.run([sys.executable, '-m', 'venv', dir_venv], check=True)
    subprocess.run([py, '-m', 'pip', 'install'] + requirements, check=True)
    token = input('Enter your bot token or nothing to start a preconfigured bot.')
    if token == '':
        with open(file_json, 'r') as f:
            data = json.load(f)
            token = data['token']
        subprocess.run([py, file_py, token], check=True)
    else:
        print('Configuring....')
        if not os.path.exists(dir_json):
            os.mkdir(dir_json)
        if not os.path.exists(dir_log):
            os.mkdir(dir_log)
        if not os.path.exists(dir_txt):
            os.mkdir(dir_txt)
        with open(file_json, 'w') as f:
            json.dump({
                'token': token,
                'prefix': "!",
                'os': system
            }, f, indent=4)
        subprocess.run([py, file_py, token], check=True)
