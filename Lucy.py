import os
import subprocess
import sys
import venv

# Directory paths
base_dir = os.path.dirname(os.path.abspath(__file__))
venv_dir = os.path.join(base_dir, 'venv')

# Requirements
requirements = [
    'discord.py',
    'rdkit',
    'Pillow',
    'pyyaml',
    'requests'
]

# Create virtual environment if it does not exist
if not os.path.exists(venv_dir):
    print("Creating virtual environment...")
    venv.create(venv_dir, with_pip=True)
else:
    print("Virtual environment already exists.")

# Get the path to the Python interpreter in the virtual environment
if os.name == 'nt':
    venv_python = os.path.join(venv_dir, 'Scripts', 'python.exe')
else:
    venv_python = os.path.join(venv_dir, 'bin', 'python')

# Function to run a command and handle errors
def run_command(command):
    try:
        subprocess.check_call(command)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running command: {' '.join(command)}")
        print(e)
        sys.exit(1)

# Install packages
print("Installing required packages...")
run_command([venv_python, '-m', 'pip', 'install', '--upgrade', 'pip'])
run_command([venv_python, '-m', 'pip', 'install'] + requirements)

# Run the main bot script using the Python interpreter from the virtual environment
main_script = os.path.join(base_dir, 'Main.py')
print("Running the bot script...")

run_command([venv_python, main_script])
