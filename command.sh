#!/bin/bash

# Directory of the script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Name of the virtual environment
VENV_DIR="$SCRIPT_DIR/venv"

# Required packages
REQUIRED_PACKAGES=("discord.py" "rdkit" "requests" "pillow" "pubchempy" "pyyaml")

# Check if virtual environment exists, if not create it
if [ ! -d "$VENV_DIR" ]; then
  echo "Creating virtual environment..."
  python3 -m venv "$VENV_DIR"
fi

# Activate the virtual environment
source "$VENV_DIR/bin/activate"

# Install required packages if not already installed
for package in "${REQUIRED_PACKAGES[@]}"; do
  if ! pip show "$package" > /dev/null 2>&1; then
    echo "Installing $package..."
    pip install "$package"
  else
    echo "$package is already installed."
  fi
done

# Execute main.py
python "$SCRIPT_DIR/Main.py"

# Deactivate the virtual environment
deactivate

echo "Setup complete. Virtual environment is ready and required packages are installed."

