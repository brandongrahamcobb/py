""" helpers.py
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

from discord.ext import commands
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw

import discord
import json
import os
import rdkit

up_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
up_up_dir = os.path.join(up_dir, '..')

def monograph(mol: rdkit.Chem.rdchem.Mol) -> discord.File:
    img = Draw.MolToImage(mol, size=(512, 512))
    img_bytes = io.BytesIO()
    img.save(img_bytes, format='PNG')
    img_bytes.seek(0)
    file = discord.File(fp=img_bytes, filename=f'Molecule.png')
    return file

def get_version():
    if not os.path.exists(os.path.join(up_up_dir, 'version.txt')):
        with open(os.path.join(up_up_dir, 'version.txt'), 'w') as f:
            f.write('1.0.0')
    with open(os.path.join(up_up_dir, 'version.txt'), 'r') as f:
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
    with open(os.path.join(up_up_dir, 'version.txt'), 'w') as f:
        f.write(version)
    return version

def load_config():
    config_path = os.path.join(up_up_dir, 'config.json')
    if os.path.exists(config_path):
        with open(config_path, 'r') as f:
            return json.load(f)
    else:
        raise FileNotFoundError("Configuration file not found.")
