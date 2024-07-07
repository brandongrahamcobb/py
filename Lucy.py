import os
import subprocess
import sys
import platform

# Directory paths
base_dir = os.path.dirname(os.path.abspath(__file__))
venv_dir = os.path.join(base_dir, 'venv')
log_dir = os.path.join(base_dir, '..', 'log')
txt_dir = os.path.join(base_dir, '..', 'txt')

# Requirements
requirements = [
    'discord.py',
    'rdkit',
    'pillow',
    'pyyaml',
    'requests'
]

# Function to run a command and handle errors
def run_command(command):
    try:
        subprocess.check_call(command)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running command: {' '.join(command)}")
        print(e)
        sys.exit(1)

# Determine the operating system
os_type = platform.system()

if os_type == 'Windows':
    # Install packages globally on Windows
    print("Detected Windows OS. Installing required packages globally...")
    for package in requirements:
        run_command([sys.executable, '-m', 'pip', 'install', '--upgrade', package])
    
    # Run the main bot script using the system Python interpreter
    main_script = os.path.join(base_dir, 'Main.py')
    print("Running the bot script...")
    run_command([sys.executable, main_script])

    if not os.path.exists(log_dir):
        print("Creating log folder...")
        os.makedirs(log_dir)
    else:
        print("Log folder already exists.")

    if not os.path.exists(txt_dir):
        print("Creating txt folder...")
        os.makedirs(txt_dir)
    else:
        print("Text folder already exists.")
else:
    # Create virtual environment if it does not exist
    if not os.path.exists(venv_dir):
        print("Creating virtual environment...")
        venv.create(venv_dir, with_pip=True)
    else:
        print("Virtual environment already exists.")

    if not os.path.exists(log_dir):
        print("Creating log folder...")
        os.makedirs(log_dir)
    else:
        print("Log folder already exists.")

    if not os.path.exists(txt_dir):
        print("Creating txt folder...")
        os.makedirs(txt_dir)
    else:
        print("Text folder already exists.")

    # Get the path to the Python interpreter in the virtual environment
    venv_python = os.path.join(venv_dir, 'bin', 'python') if os_type == 'Linux' or os_type == 'Darwin' else os.path.join(venv_dir, 'Scripts', 'python.exe')

    # Install packages
    print("Installing required packages...")
    run_command([venv_python, '-m', 'pip', 'install', '--upgrade', 'pip'])
    run_command([venv_python, '-m', 'pip', 'install'] + requirements)

    # Run the main bot script using the Python interpreter from the virtual environment
    main_script = os.path.join(base_dir, 'Main.py')
    print("Running the bot script...")
    run_command([venv_python, main_script])
