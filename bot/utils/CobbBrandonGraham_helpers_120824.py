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

from googleapiclient.discovery import build
def google_search(search_term, api_key, cse_id, **kwargs):
    service = build('customsearch', 'v1', developerKey=api_key)
    res = service.cse().list(q=search_term, cx=cse_id, **kwargs).execute()
    return res['items']
def get_cse_id(option):
    cse_ids = {
        'turkey': '25fb95395119b40f0',
        'web': 'e1753e22095c74f10',
        'msds': '4231a573aa34241e4',
        'fda': '72a76e05b0ba044d0'
    }
    return cse_ids.get(option, cse_ids['web'])
def search(query, option):
    config = load_config()
    google_json_key = config['api_keys']['api_key_1']
    cse_id = get_cse_id(option)
    results = google_search(query, google_json_key, cse_id, num=3)
    if results:
       return results
    return None

import emoji as emoji_lib
import unicodedata
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

import pubchempy as pcp
def get_sds(query: str):
    compound = get_molecule_info(query)
    compound_info = compound.to_dict(properties=['iupac_name', 'molecular_formula', 'molecular_weight', 'synonyms'])
    description = compound_info.get('synonyms', ['No description available'])[0]
    pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{compound.cid}"
    fda_guide_url = search(query, 'fda')
    msds_message = (
      f"**Molecule Name:** {compound_info.get('iupac_name', 'N/A')}\n"
      f"**Molecular Formula:** {compound_info.get('molecular_formula', 'N/A')}\n"
      f"**Molecular Weight:** {compound_info.get('molecular_weight', 'N/A')} g/mol\n"
      f"**Description:** {description}\n"
      f"**PubChem URL:** [View PubChem]({pubchem_url})\n"
      f"**FDA Prescriber Guide:** {fda_guide_url[0]['link'] if fda_guide_url[0]['link'] else 'No prescriber guide found.'}"
    )
    return msds_message

def get_molecule_info(query: str):
    compounds = pcp.get_compounds(query, 'name')
    return compounds[0]

from PIL import Image, ImageDraw, ImageFont

def calculate_max_font_size(text, max_width, min_size=10, max_size=30):
    font_size = max_size
    temp_image = Image.new("RGB", (0, 0))
    draw = ImageDraw.Draw(temp_image)
    while font_size >= min_size:
        font = ImageFont.load_default(font_size)
        text_bbox = draw.textbbox((0, 0), text, font=font)
        text_width = text_bbox[2] - text_bbox[0]
        if text_width <= max_width:
            return font_size
        font_size -= 1
    return min_size
