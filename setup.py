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
import shutil

base_dir = os.path.abspath(os.path.dirname(__file__))
py_dir = os.path.join(base_dir, 'py')

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

    venv_dir = os.path.join(base_dir, 'venv')
    
    if not os.path.exists(venv_dir):
        # Create virtual environment
        subprocess.run([sys.executable, '-m', 'venv', venv_dir], check=True)
    
    requirements = ["discord.py", "emoji", "pubchempy", "rdkit", "pillow", "requests"]
    python_exec = os.path.join(venv_dir, 'bin', 'python') if system == 'Linux' else os.path.join(venv_dir, 'Scripts', 'python')
    
    # Update pip and install requirements
    subprocess.run([python_exec, '-m', 'pip', 'install', '--upgrade', 'pip'], check=True)
    subprocess.run([python_exec, '-m', 'pip', 'install', '--upgrade'] + requirements, check=True)

    # Create necessary directories at the base directory level
    os.makedirs(os.path.join(base_dir, 'json'), exist_ok=True)
    os.makedirs(os.path.join(base_dir, 'log'), exist_ok=True)
    os.makedirs(os.path.join(base_dir, 'txt'), exist_ok=True)

    # Move this setup script, main.py, my_cog.py, and game_cog.py to py directory if not already
    setup_script_path = os.path.join(base_dir, 'setup.py')
    main_script_path = os.path.join(base_dir, 'main.py')
    my_cog_script_path = os.path.join(base_dir, 'my_cog.py')
    game_cog_script_path = os.path.join(base_dir, 'game_cog.py')
    license_path = os.path.join(base_dir, 'LICENSE')
    readme_path = os.path.join(base_dir, 'README.md')
    
    if not os.path.exists(py_dir):
        os.makedirs(py_dir)

    if os.path.exists(main_script_path):
        shutil.move(main_script_path, os.path.join(py_dir, 'main.py'))

    if os.path.exists(my_cog_script_path):
        shutil.move(my_cog_script_path, os.path.join(py_dir, 'my_cog.py'))

    if os.path.exists(game_cog_script_path):
        shutil.move(game_cog_script_path, os.path.join(py_dir, 'game_cog.py'))

    if os.path.exists(license_path):
        shutil.move(license_path, os.path.join(py_dir, 'LICENSE'))

    if os.path.exists(readme_path):
        shutil.move(readme_path, os.path.join(py_dir, 'README.md'))
    
    new_setup_script_path = os.path.join(py_dir, 'setup.py')
    if os.path.abspath(__file__) != new_setup_script_path:
        shutil.move(setup_script_path, new_setup_script_path)
        print("Moved setup script to 'py' directory.")
        return  # Exit to prevent running the script twice

def set_token():
    """Prompt user to input bot token if not present."""
    config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'json', 'config.json')
    
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
    if os.path.basename(os.path.dirname(__file__)) == 'py':
        print("Setup script is inside 'py' directory. Doing nothing.")
    else:
        setup_environment()
        set_token()
        print(f"Setup complete. New version: {get_version()}")
