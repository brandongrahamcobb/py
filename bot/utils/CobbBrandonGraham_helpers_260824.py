''' helpers.py
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
'''

from bot.main import Lucy
from datetime import datetime
from discord.ext import commands
from functools import wraps
from googleapiclient.discovery import build
from os.path import isfile
from rdkit import Chem

import asyncio
import discord
import emoji as emoji_lib
import io
import itertools
import os
import pubchempy as pcp
import random
import requests
import unicodedata
import json
import yaml
import wikipediaapi

def get_emoji(emoji_character):
    if emoji_character is None:
        return 'Please provide a Unicode emoji character.'
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
        return unicode_info
    else:
        return 'Unicode emoji not found. Make sure it is a valid Unicode emoji character.'

def get_images(query, num_results=30):
    service = build('customsearch', 'v1', developerKey=(load_config()['api_keys']['api_key_1']))
    image_urls = []
    start_index = 1
    while len(image_urls) < num_results:
        response = service.cse().list(
            q=query,
            cx='e1753e22095c74f10',
            searchType='image',
            num=10,
            start=start_index
        ).execute()
        current_images = [item['link'] for item in response.get('items', [])]
        if not current_images:
            break  # Exit loop if no more images are found
        image_urls.extend(current_images)
        start_index += 10  # Increment to get the next set of results
    if not image_urls:
        raise ValueError('No images found.')
    return random.choice(image_urls)

def get_mol(arg):
    compounds = pcp.get_compounds(arg, 'name') #, record_type='3d')  # '3d' to get the most relevant record
    compound_data = compounds[0].to_dict(properties=['isomeric_smiles'])
    mol = Chem.MolFromSmiles(compound_data['isomeric_smiles'])
    if mol is None:
        raise ValueError('Invalid SMILES string')
    return mol


def get_sds(query: str):
    # Convert query to lowercase for easier comparison
    query_lower = query.lower()

    # Custom handling for the "omnibol" combination of DHEA and pregnenolone
    if query_lower == "omnibol":
        description = "Omnibol represents the combination of DHEA and Pregnenolone."
        dhea_info = pcp.get_compounds("dhea", 'name')[0].to_dict(properties=['iupac_name', 'molecular_formula', 'molecular_weight'])
        pregnenolone_info = pcp.get_compounds("pregnenolone", 'name')[0].to_dict(properties=['iupac_name', 'molecular_formula', 'molecular_weight'])

        embed = discord.Embed(title="Omnibol (DHEA + Pregnenolone)", description=description, color=0x00ff99)
        embed.add_field(name="DHEA Name", value=dhea_info.get('iupac_name', 'N/A'), inline=True)
        embed.add_field(name="DHEA Molecular Formula", value=dhea_info.get('molecular_formula', 'N/A'), inline=True)
        embed.add_field(name="DHEA Molecular Weight", value=f"{dhea_info.get('molecular_weight', 'N/A')} g/mol", inline=True)

        embed.add_field(name="Pregnenolone Name", value=pregnenolone_info.get('iupac_name', 'N/A'), inline=True)
        embed.add_field(name="Pregnenolone Molecular Formula", value=pregnenolone_info.get('molecular_formula', 'N/A'), inline=True)
        embed.add_field(name="Pregnenolone Molecular Weight", value=f"{pregnenolone_info.get('molecular_weight', 'N/A')} g/mol", inline=True)

        embed.add_field(name="Your Input", value=query, inline=True)
        return embed

    # Custom handling for DHEA and Pregnenolone individually
    if query_lower in ["dhea", "pregnenolone"]:
        compound_info = pcp.get_compounds(query_lower, 'name')[0].to_dict(properties=['iupac_name', 'molecular_formula', 'molecular_weight', 'synonyms'])
        description = f"{query.capitalize()} is a steroid precursor with hormone-regulating functions."
        embed = discord.Embed(title=f"Custom {query.capitalize()} Information", description=description, color=0x00ff99)
        embed.add_field(name="Molecule Name", value=compound_info.get('iupac_name', 'N/A'), inline=False)
        embed.add_field(name="Your Input", value=query, inline=True)
        embed.add_field(name="Additional Info", value="For more details, consult medical references.", inline=True)
        return embed

    # Standard handling for other compounds
    try:
        # Fetch the compound based on the query
        compound = pcp.get_compounds(query, 'name')[0]
        compound_info = compound.to_dict(properties=['iupac_name', 'molecular_formula', 'molecular_weight', 'synonyms'])
        description = compound_info.get('synonyms', ['No description available'])[0]
        pubchem_url = f'https://pubchem.ncbi.nlm.nih.gov/compound/{compound.cid}'
        pubmed_description_url = f'https://www.ncbi.nlm.nih.gov/sites/entrez?LinkName=pccompound_pubmed&db=pccompound&cmd=Link&from_uid={compound.cid}'

        # Fetch the PubChem description
        description = compound.record['description'][0] if 'description' in compound.record else description
    except (KeyError, IndexError):
        description = 'No detailed description available.'

    embed = discord.Embed(title=description, color=0x00ff99)
    embed.add_field(name="Molecule Name", value=compound_info.get('iupac_name', 'N/A'), inline=False)
    embed.add_field(name="Your Input", value=query, inline=True)
    embed.add_field(name="PubChem URL", value=f"[View PubChem]({pubchem_url})" if pubchem_url else "No PubChem URL found.", inline=True)
    embed.add_field(name="PubMed URL", value=f"[View PubMed]({pubmed_description_url})" if pubmed_description_url else "No PubMed URL found.", inline=True)

    return embed

def is_off_peak_hours():
    '''
    Determine if the current time is during off-peak hours.
    Off-peak hours: Weekdays from 9:00 PM to 5:00 AM ET and weekends.
    '''
    now = datetime.now()
    if now.weekday() >= 5:  # Saturday (5) or Sunday (6)
        return True
    if now.hour > 5 or now.hour <= 21:  # Before 5 AM or after 9 PM
        return True
    return False

def load_config():
    home_dir = os.path.expanduser('~')
    config_path = os.path.join(home_dir, '.config', 'lucy', 'config.json')
    if not os.path.exists(config_path):
        raise FileNotFoundError('Configuration file not found.')
    with open(config_path, 'r') as f:
        return json.load(f)

def read_config():
    home_dir = os.path.expanduser('~')
    config_path = os.path.join(home_dir, '.config', 'lucy', 'config.json')
    if not os.path.exists(config_path):
        raise FileNotFoundError('Configuration file not found.')
    with open(config_path, 'r') as f:
        return json.load(f)

def read_users():
    home_dir = os.path.expanduser('~')
    users_path = os.path.join(home_dir, '.config', 'lucy', 'users.yaml')
    if not os.path.exists(config_path):
        raise FileNotFoundError('User file not found.')
    with open(users_path, 'r') as f:
        return yaml.safe_load(f)

def write_users(data):
    home_dir = os.path.expanduser('~')
    users_path = os.path.join(home_dir, '.config', 'lucy', 'users.yaml')
    if not os.path.exists(user_path):
        raise FileNotFoundError('User file not found.')
    with open(users_path, 'w') as f:
        return yaml.dump(data, f)

def off_peak_hours_required():
    '''
    Decorator to restrict command execution to off-peak hours.
    '''
    def decorator(func):
        @wraps(func)
        async def wrapper(*args, **kwargs):
            self = args[0]  # The first argument is typically `self` or `ctx`
            ctx = args[1]  # The first argument is typically `self` or `ctx`
            if is_off_peak_hours():
                return await func(*args, **kwargs)
            else:
                await ctx.send('This operation is restricted to off-peak hours (9:00 PM - 5:00 AM ET on weekdays, or anytime on weekends). Please try again later.')
        return wrapper
    return decorator

def unique_pairs(strings_list):
    pairs = list(itertools.combinations(strings_list, 2))
    sorted_pairs = [sorted(list(pair)) for pair in pairs]
    sorted_pairs_overall = sorted(sorted_pairs)
    return sorted_pairs_overall
