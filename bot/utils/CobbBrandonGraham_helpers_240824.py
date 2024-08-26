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

from bot.main import Lucy
from datetime import datetime
from os.path import isfile
from discord.ext import commands
from functools import wraps
from rdkit import Chem
import asyncio
import discord
import emoji as emoji_lib
import io
import itertools
import json
import os
import pubchempy as pcp
import random
import requests
import unicodedata
from rdkit import Chem

from googleapiclient.discovery import build

NCBI_REQUEST_DELAY = 0.34  # Delay to ensure no more than 3 requests per second
tasks = {}  # Dictionary to store tasks keyed by (command name, user ID)

def unique_pairs(strings_list):
    pairs = list(itertools.combinations(strings_list, 2))
    sorted_pairs = [sorted(list(pair)) for pair in pairs]
    sorted_pairs_overall = sorted(sorted_pairs)
#    processed_pairs = [[repr(s)[1:-1] for s in pair] for pair in sorted_pairs_overall]
 #   return processed_pairs
    return sorted_pairs_overall

async def delayed_command(ctx, delay_seconds):
    task_key = (ctx.command.name, ctx.author.id)
    if task_key in tasks and tasks[task_key]:
        tasks[task_key].cancel()
    async def delayed_action():
        await asyncio.sleep(delay_seconds)
        def is_bot_message(message):
            return message.author == self.bot.user
        deleted = await self.purge_messages(ctx, limit, is_bot_message)
    tasks[task_key] = asyncio.create_task(delayed_action())

def fetch_ncbi_organism(tax_id):
    """
    Fetch detailed information about an organism from NCBI Taxonomy database.
    """
#    https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
    summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    summary_params = {
        "db": "taxonomy",
        "id": tax_id,
        "retmode": "json"
    }
    response = requests.get(summary_url, params=summary_params)
    data = response.json()
    if response.status_code == 200 and "result" in data and tax_id in data["result"]:
        organism_info = data["result"][tax_id]
        return {
            "tax_id": tax_id,
            "scientific_name": organism_info.get("scientificname", "N/A"),
            "rank": organism_info.get("rank", "N/A"),
            "division": organism_info.get("division", "N/A"),
            "genetic_code": organism_info.get("genetic_code", "N/A"),
            "lineage": organism_info.get("lineage", "N/A"),
            "description": organism_info.get("description", "N/A")
        }
    else:
        return {
            "tax_id": tax_id,
            "scientific_name": "N/A",
            "rank": "N/A",
            "division": "N/A",
            "genetic_code": "N/A",
            "lineage": "N/A",
            "description": "N/A"
        }

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

def get_sds(query: str):
    compound = pcp.get_compounds(query, 'name')[0]
    compound_info = compound.to_dict(properties=['iupac_name', 'molecular_formula', 'molecular_weight', 'synonyms'])
    description = compound_info.get('synonyms', ['No description available'])[0]
    pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{compound.cid}"
#    fda_guide_url = search(query, 'fda')
    msds_message = (
      f"**Molecule Name:** {compound_info.get('iupac_name', 'N/A')}\n"
      f"**Molecular Formula:** {compound_info.get('molecular_formula', 'N/A')}\n"
      f"**Molecular Weight:** {compound_info.get('molecular_weight', 'N/A')} g/mol\n"
      f"**Description:** {description}\n"
      f"**PubChem URL:** [View PubChem]{pubchem_url if pubchem_url else 'No PubChem url found.'}"
 #     f"**FDA:** {fda_guide_url[0]['link'] if fda_guide_url[0]['link'] else 'No prescriber guide found.'}"
    )
    return msds_message

def is_off_peak_hours():
    """
    Determine if the current time is during off-peak hours.
    Off-peak hours: Weekdays from 9:00 PM to 5:00 AM ET and weekends.
    """
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
        raise FileNotFoundError("Configuration file not found.")
    with open(config_path, 'r') as f:
        return json.load(f)

def off_peak_hours_required():
    """
    Decorator to restrict command execution to off-peak hours.
    """
    def decorator(func):
        @wraps(func)
        async def wrapper(*args, **kwargs):
            self = args[0]  # The first argument is typically `self` or `ctx`
            ctx = args[1]  # The first argument is typically `self` or `ctx`
            if is_off_peak_hours():
                return await func(*args, **kwargs)
            else:
                await ctx.send("This operation is restricted to off-peak hours (9:00 PM - 5:00 AM ET on weekdays, or anytime on weekends). Please try again later.")
        return wrapper
    return decorator

def search_taxonomy(term, max_results=10):
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    search_params = {
        "db": "taxonomy",
        "term": term,
        "retmode": "json",
        "retmax": max_results
    }
    search_response = requests.get(search_url, params=search_params)
    search_data = search_response.json()
    if not search_data["esearchresult"]["idlist"]:
        return []
    tax_ids = search_data["esearchresult"]["idlist"]
    return [{"tax_id": tax_id} for tax_id in tax_ids]

def validate_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES")
        return True
    except Exception as e:
        print(f"Error validating SMILES '{smiles}': {e}")
        return False

def get_mol(arg):
    compounds = pcp.get_compounds(arg, 'name') #, record_type='3d')  # '3d' to get the most relevant record
    compound_data = compounds[0].to_dict(properties=['isomeric_smiles'])
    mol = Chem.MolFromSmiles(compound_data['isomeric_smiles'])
    if mol is None:
        raise ValueError("Invalid SMILES string")
    return mol

def get_images(query, num_results=30):
    service = build("customsearch", "v1", developerKey=(load_config()['api_keys']['api_key_1']))
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
        raise ValueError("No images found.")
    return random.choice(image_urls)
