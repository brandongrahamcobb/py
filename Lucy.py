# File: Lucy.py
# Description: The purpose of this file is to serve as an entry-point for Lucy, the discord.py bot. The program has a first run function and a default run function.
# Author: spawd
# Author: Lily
# Date: July 6, 2024
# Usage: python Lucy.py
# Notes: <Any additional notes or information about the file>
# Dependencies: Chemistry.py Listener.py Main.py Super.py
'''
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
import venv

requirements = ["discord.py", "pubchempy", "rdkit", "pillow", "requests"]

dir = os.path.dirname(os.path.abspath(__file__))
dir_json = os.path.join(dir, 'json')
dir_log = os.path.join(dir, 'log')
dir_py = os.path.join(dir, 'py')
dir_txt = os.path.join(dir, 'txt')
dir_venv = os.path.join(dir, 'venv')
file_json = os.path.join(dir_json, 'config.json')
file_py_0 = os.path.join(dir, 'Main.py')
file_py_1 = os.path.join(dir, 'Chemistry.py')
file_py_2 = os.path.join(dir, 'Listener.py')
file_py_3 = os.path.join(dir, 'Super.py')
file_py_4 = os.path.join(dir_py, 'Main.py')
file_py_5 = os.path.join(dir_py, 'Chemistry.py')
file_py_6 = os.path.join(dir_py, 'Listener.py')
file_py_7 = os.path.join(dir_py, 'Super.py')
system = 'Linux' if os.name != 'nt' else 'Windows'
py = os.path.join(dir_venv, 'bin', 'python') if system == 'Linux' else sys.executable


def on_ready():
    print(f'System: {system}')
    print(f'Location: {dir}')
    print(f'Python Location: {sys.executable}')
    print(f'Python Version: {sys.version}')

def configuration():
    token = input('Enter your bot token.')
    if not os.path.isdir(dir_json):
        os.mkdir(dir_json)
    if not os.path.isdir(dir_py):
        os.mkdir(dir_py)
        os.replace(file_py_0, file_py_4)
        os.replace(file_py_1, file_py_5)
        os.replace(file_py_2, file_py_6)
        os.replace(file_py_3, file_py_7)
    if not os.path.isdir(dir_log):
        os.mkdir(dir_log)
    if not os.path.isdir(dir_txt):
        os.mkdir(dir_txt)
    if system == 'Linux':
        if not os.path.isdir(dir_venv):
            venv.create(dir_venv, with_pip=True)
        subprocess.run([py, '-m', 'venv', 'activate'], check=True)
    subprocess.run([py, '-m', 'pip', 'install'] + requirements, check=True)
    if os.path.isfile(file_json):
        os.remove(file_json)
    with open(file_json, 'w') as f:
        json.dump({
            'token': token,
            'prefix': "!",
            'os': system
        }, f, indent=4)
    print('Done configuring...')
    with open(file_json, 'r') as f:
        data = json.load(f)
        token = data['token']
    subprocess.run(py + ' ' + file_py_4 + ' ' + token, check=True, shell=True)

if __name__ == '__main__':
    try:
        on_ready()
        configuration()
    except Exception as e:
        None
