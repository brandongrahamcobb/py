""" setup.py
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
"""

import os
import subprocess
import sys
import json

def is_root():
    return os.geteuid() == 0  # Check if running as root

def setup_environment():
    """Setup the environment, including virtual environment and requirements."""
    system = 'Linux' if os.name != 'nt' else 'Windows'
    
    if system == 'Linux' and is_root():
        print('Running as root, updating system...')
        subprocess.run(['pacman', '-Syu'], check=True)
    elif system == 'Linux':
        print('Not running as root, skipping system update.')

    venv_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'venv')
    if not os.path.exists(venv_dir):
        # Create virtual environment
        subprocess.run([sys.executable, '-m', 'venv', venv_dir], check=True)
    
    requirements = ["asyncpraw", "discord.py", "emoji", "pubchempy", "rdkit", "pillow", "requests"]
    python_exec = os.path.join(venv_dir, 'bin', 'python') if system == 'Linux' else sys.executable
    subprocess.run([python_exec, '-m', 'pip', 'install'] + requirements, check=True)

    # Create necessary directories
    os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'json'), exist_ok=True)
    os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'log'), exist_ok=True)
    os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'txt'), exist_ok=True)

def set_token():
    """Prompt user to input bot token if not present."""
    config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'json', 'config.json')
    
    if not os.path.exists(config_path):
        token = input("Enter your Discord bot token: ")
        os.makedirs(os.path.dirname(config_path), exist_ok=True)
        with open(config_path, 'w') as f:
            json.dump({"token": token}, f, indent=4)
        print("Token has been saved to config.json.")
    else:
        with open(config_path, 'r') as f:
            config = json.load(f)
        if "token" not in config:
            token = input("Enter your Discord bot token: ")
            config["token"] = token
            with open(config_path, 'w') as f:
                json.dump(config, f, indent=4)
            print("Token has been updated in config.json.")
        else:
            print("Token already set in config.json.")

def get_version():
    version_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'txt', 'version.txt')
    if not os.path.exists(version_file):
        with open(version_file, 'w') as f:
            f.write('1.0.0')
    with open(version_file, 'r') as f:
        version = f.read().strip()
    major, minor, patch = map(int, version.split('.'))
    patch += 1
    if patch >= 10:
        patch = 0
        minor += 1
    if minor >= 10:
        minor = 0
        major += 1
    version = f"{major}.{minor}.{patch}"
    with open(version_file, 'w') as f:
        f.write(version)
    return version

if __name__ == '__main__':
    setup_environment()
    set_token()
    print(f"Setup complete. New version: {get_version()}")
