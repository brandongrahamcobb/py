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
file_py = os.path.join(dir, 'Main.py')
system = 'Linux' if os.name != 'nt' else 'Windows'
py = os.path.join(dir_venv, 'bin', 'python') if system == 'Linux' else sys.executable

if __name__ == '__main__':
    print(f'System: {system}')
    print(f'Location: {dir}')
    print(f'Python Location: {sys.executable}')
    print(f'Python Version: {sys.version}')
    if not os.path.isdir(dir_json) and not os.path.isdir(dir_log) and not os.path.isdir(dir_txt) and not os.path.isdir(dir_venv):
        print(f'Configuring....')
        os.mkdir(dir_json)
        os.mkdir(dir_log)
        os.mkdir(dir_txt)
        token = input('Enter your bot token.')
        with open(file_json, 'w') as f:
            json.dump({
                'token': token,
                'prefix': "!",
                'os': system
            }, f, indent=4)
        subprocess.run([sys.executable, '-m', 'venv', dir_venv], check=True)
    print('Updating packages...')
    subprocess.run([py, '-m', 'pip', 'install'] + requirements, check=True)
    print('Done updating...')
    with open(file_json, 'r') as f:
        data = json.load(f)
        token = data['token']
    subprocess.run(py + ' ' + file_py + ' ' + token, check=True, shell=True)
