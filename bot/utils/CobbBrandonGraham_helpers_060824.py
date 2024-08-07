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

import logging
import logging.handlers
def setup_logging():
    global logger
    log_file = 'discord.log'
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.addHandler(file_handler)
    if not os.path.exists(log_file):
        open(log_file, 'a').close()

import json
import os
def load_config():
    config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config.json')
    if not os.path.exists(config_path):
        raise FileNotFoundError("Configuration file not found.")
    with open(config_path, 'r') as f:
        return json.load(f)
def get_version():
    version_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'version.txt')
    try:
        with open(version_file, 'r') as f:
            version = f.read().strip()
    except FileNotFoundError:
        version = None
    return version

import requests
from bs4 import BeautifulSoup
def fetch_and_parse(url):
    try:
        response = requests.get(url)
        response.raise_for_status()
        soup = BeautifulSoup(response.content, 'html.parser')
        return soup.get_text()
    except requests.RequestException as e:
        print(f"Request failed: {e}")
        return None

import hashlib
def hash_content(content):
    return hashlib.sha256(content.encode()).hexdigest()
def same_hashes(url1, url2):
    try:
        response1 = requests.get(url1)
        response2 = requests.get(url2)
        if response1.status_code == 200 and response2.status_code == 200:
            hash1 = hash_content(response1.text)
            hash2 = hash_content(response2.text)
            if hash1 == hash2:
                return True
            else:
                return False
    except requests.RequestException as e:
        return f"An error occurred: {e}"

from googleapiclient.discovery import build
def google_search(search_term, api_key, cse_id, **kwargs):
    service = build('customsearch', 'v1', developerKey=api_key)
    res = service.cse().list(q=search_term, cx=cse_id, **kwargs).execute()
    return res['items']

import emoji as emoji_lib
import unicodedata
async def get_emoji_info(ctx, emoji_character):
    if emoji_character is None:
        await ctx.send('Please provide a Unicode emoji character.')
        return
    unicode_name = emoji_lib.demojize(emoji_character)
    if unicode_name.startswith(':') and unicode_name.endswith(':'):
        unicode_name = unicode_name[1:-1]
        code_points = ' '.join(f'U+{ord(c):04X}' for c in emoji_character)
        description = unicodedata.name(emoji_character, 'No description available')
        unicode_info = (
            f'**Unicode Emoji Character:** {emoji_character}\n'
            f'**Unicode Emoji Name:** {unicode_name}\n'
            f'**Code Points:** {code_points}\n'
            f'**Description:** {description}'
        )
        await ctx.send(unicode_info)
    else:
        await ctx.send('Unicode emoji not found. Make sure it is a valid Unicode emoji character.')

import pubchempy as pcp
async def get_msds_info(ctx, molecule_name: str):
    compounds = pcp.get_compounds(molecule_name, 'name')
    if not compounds:
        await ctx.send('Molecule not found.')
        return
    compound = compounds[0]
    compound_info = compound.to_dict(properties=['iupac_name', 'molecular_formula', 'molecular_weight', 'synonyms'])
    description = compound_info.get('synonyms', ['No description available'])[0]
    msds_url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{compound.cid}"
#    fda_guide_url = self.get_fda_prescriber_guide(molecule_name)
    msds_message = (
      f"**Molecule Name:** {compound_info.get('iupac_name', 'N/A')}\n"
      f"**Molecular Formula:** {compound_info.get('molecular_formula', 'N/A')}\n"
      f"**Molecular Weight:** {compound_info.get('molecular_weight', 'N/A')} g/mol\n"
      f"**Description:** {description}\n"
      f"**MSDS URL:** [View MSDS]({msds_url})\n"
 #     f"**FDA Prescriber Guide:** {fda_guide_url if fda_guide_url else 'No prescriber guide found.'}"
    )
    await ctx.send(msds_message)
