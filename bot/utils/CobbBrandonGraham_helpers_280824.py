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
from googleapiclient.discovery import build
from io import BytesIO
from os import makedirs
from os.path import abspath, dirname, exists, expanduser, isfile, join
from PIL import Image, ImageFont, ImageDraw
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, DataStructs, Draw
from rdkit.Chem.Draw import rdMolDraw2D, SimilarityMaps
from rdkit.DataStructs import FingerprintSimilarity, TanimotoSimilarity
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from typing import List
from webdriver_manager.chrome import ChromeDriverManager
from typing import List, Optional, Dict, Any

import asyncio
import bot.utils.helpers as help
import colorsys
import datetime as dt
import discord
import emoji as emoji_lib
import itertools
import logging
import logging.handlers
import os
import pubchempy as pcp
import math
import openai
import random
import requests
import unicodedata
import yaml

current_date = dt.datetime.now().strftime('%d%m%y')

home = expanduser('~')


path_base = join(home, 'Documents', 'src', 'lucy')
path_ai_cog = join(path_base, 'cogs', 'ai_cog.py')
path_config_yaml = join(home, '.config', 'lucy', 'config.yaml')
path_log = join(home, '.log', 'discord.log')
path_users_yaml = join(home, '.config', 'lucy', 'users.yaml')

def add_watermark(image: BytesIO, watermark_text: str = '~LUCY~') -> BytesIO:
    RGB_image = Image.open(image)
    RGBA_image = RGB_image.convert('RGBA')
    draw = ImageDraw.Draw(RGBA_image)
    width, height = RGBA_image.size
    diagonal = math.sqrt(width**2 + height**2)
    font_size = int(diagonal / 15)
    try:
        font = ImageFont.truetype('resources/Roboto-Regular.ttf', font_size)  # Replace with the path to your font file
    except IOError:
        font = ImageFont.load_default()
    while True:
        bbox = draw.textbbox((0, 0), watermark_text, font=font)
        text_width = bbox[2] - bbox[0]
        if text_width <= 512 or font_size <= 1:
            break
        font_size -= 1
        font = ImageFont.truetype('resources/Roboto-Regular.ttf', font_size)  # Replace with the path to your font file
    text_height = bbox[3] - bbox[1]
    text_x = (width - text_width) / 2
    text_y = height - (2 * text_height)
    watermark_image = Image.new('RGBA', RGBA_image.size, (0, 0, 0, 0))
    watermark_draw = ImageDraw.Draw(watermark_image)
    watermark_draw.text((text_x, text_y), watermark_text, font=font, fill=(255, 255, 255, 64))
    mask = watermark_image.split()[3]
    RGBA_image.paste(watermark_image, (0, 0), mask)
    output = BytesIO()
    RGBA_image.save(output, format= 'PNG')
    output.seek(0)
    return output

def adjust_hue_and_saturation(image, hue_shift, saturation_shift) -> BytesIO:
    image = image.convert('RGB')
    pixels = list(image.getdata())
    adjusted_pixels = []
    for r, g, b in pixels:
        h, s, v = colorsys.rgb_to_hsv(r / 255.0, g / 255.0, b / 255.0)
        h = (h + hue_shift / 360.0) % 1.0
        s = min(max(s + saturation_shift / 100.0, 0), 1)
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        adjusted_pixels.append((int(r * 255), int(g * 255), int(b * 255)))
    new_image = Image.new('RGB', image.size)
    new_image.putdata(adjusted_pixels)
    output = BytesIO()
    new_image.save(output, format= 'PNG')
    output.seek(0)
    return output

async def check_for_updates(bot):
    openai.api_key = bot.config['api_keys']['api_key_1']
    while True:
        await asyncio.sleep(360)  # Wait for 6 minutes
        prompt = 'Please provide the updated code for the AICog class. Ensure the code is within the scope of a cog and does not affect other parts of the bot.'
        response = openai.Completion.create(
            engine='text-davinci-003',  # Use GPT-4 if available
            prompt=prompt,
            max_tokens=500,  # Adjust this as needed
            n=1,
            stop=None,
            temperature=0.7,
        )
        generated_code = response.choices[0].text.strip()
        if validate_generated_code(generated_code):
            with open(path_ai_cog, 'r') as f:
                current_code = f.read()
            if generated_code != current_code:
                with open(path_ai_cog, 'w') as f:
                    f.write(generated_code)
                try:
                    await bot.reload_extension('bot.cogs.ai_cog')
                    print('AICog has been updated and reloaded successfully.')
                except Exception as e:
                    print(f'Failed to reload AICog: {e}')
            else:
                print('No updates needed.')
        else:
            print('Generated code did not pass validation and was not applied.')

def combine(bytes1: BytesIO, bytes2: BytesIO) -> BytesIO:
    img1 = Image.open(bytes1)
    img2 = Image.open(bytes2)
    widths, heights = zip(*(img.size for img in [img1, img2]))
    total_width = sum(widths)
    max_height = max(heights)
    combined_img = Image.new('RGB', (total_width, max_height))
    x_offset = 0
    for img in [img1, img2]:
        combined_img.paste(img, (x_offset, 0))
        x_offset += img.width
    inverted_image = invert_colors(combined_img)
    output = adjust_hue_and_saturation(inverted_image, hue_shift=-180, saturation_shift=160)
    output.seek(0)
    return output

def draw_fingerprint(pair) -> BytesIO:
    d2d = rdMolDraw2D.MolDraw2DCairo(1024, 1024)
    d2d.prepareMolsBeforeDrawing = False
    Options = d2d.drawOptions()
    Options.prepareMolsBeforeDrawing = False
    Options.includeMetadata = False
    d2d.SetDrawOptions(Options)
    mol1 = rdMolDraw2D.PrepareMolForDrawing(pair[0], kekulize=True)
    mol1.UpdatePropertyCache(False)
    mol2 = rdMolDraw2D.PrepareMolForDrawing(pair[1], kekulize=True)
    mol2.UpdatePropertyCache(False)
    fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(mol1, mol2, lambda m, i: SimilarityMaps.GetMorganFingerprint(m, i, radius=2, fpType='bv', nBits=8192), draw2d=d2d, drawingOptions=Options)
    d2d.FinishDrawing()
    drawing = d2d.GetDrawingText()
    output = BytesIO(drawing) #output = self.add_watermark(BytesIO(drawing), watermark_text=resolved_name)
    return output

def draw_watermarked_molecule(molecule) -> BytesIO:
    resolved_name = get_molecule_name(molecule)
    d2d = rdMolDraw2D.MolDraw2DCairo(1024, 1024)
    Options = d2d.drawOptions()
    d2d.SetDrawOptions(Options)
    rdMolDraw2D.SetDarkMode(Options)
    mol = rdMolDraw2D.PrepareMolForDrawing(molecule, kekulize=True)
    mol.UpdatePropertyCache(False)
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    drawing = d2d.GetDrawingText()
    output = add_watermark(BytesIO(drawing), watermark_text=resolved_name)
    return output

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
            break
        image_urls.extend(current_images)
        start_index += 10
    if not image_urls:
        raise ValueError('No images found.')
    return random.choice(image_urls)

def get_proximity(default, input) -> float:
    default_fp = AllChem.GetMorganFingerprintAsBitVect(default, 2)
    input_fp = AllChem.GetMorganFingerprintAsBitVect(input, 2)
    similarity = DataStructs.FingerprintSimilarity(default_fp, input_fp)
    return similarity

def get_mol(arg):
    compounds = pcp.get_compounds(arg, 'name')
    compound_data = compounds[0].to_dict(properties=['isomeric_smiles'])
    mol = Chem.MolFromSmiles(compound_data['isomeric_smiles'])
    if mol is None:
        raise ValueError('Invalid SMILES string')
    return mol

def get_molecule_name(molecule) -> str:
    smiles = Chem.MolToSmiles(molecule)
    if smiles:
        compounds = pcp.get_compounds(smiles, 'smiles') #, record_type='3d')
        compound_data = compounds[0].to_dict(properties=['synonyms'])
        if not compounds:
            raise ValueError('No compound found for the given SMILES string')
        return compound_data['synonyms'][0]

def get_sds(query: str):
    query_lower = query.lower()
    try:
        compound = pcp.get_compounds(query, 'name')[0]
        compound_info = compound.to_dict(properties=['iupac_name', 'molecular_formula', 'molecular_weight', 'synonyms'])
        description = compound_info.get('synonyms', ['No description available'])[0]
        pubchem_url = f'https://pubchem.ncbi.nlm.nih.gov/compound/{compound.cid}'
        pubmed_description_url = f'https://www.ncbi.nlm.nih.gov/sites/entrez?LinkName=pccompound_pubmed&db=pccompound&cmd=Link&from_uid={compound.cid}'
        description = compound.record['description'][0] if 'description' in compound.record else description
    except (KeyError, IndexError):
        description = 'No detailed description available.'
    embed = discord.Embed(title=description, color=0x00ff99)
    embed.add_field(name='Molecule Name', value=compound_info.get('iupac_name', 'N/A'), inline=False)
    embed.add_field(name='Your Input', value=query, inline=True)
    embed.add_field(name='PubChem URL', value=f'[View PubChem]({pubchem_url})' if pubchem_url else 'No PubChem URL found.', inline=True)
    embed.add_field(name='PubMed URL', value=f'[View PubMed]({pubmed_description_url})' if pubmed_description_url else 'No PubMed URL found.', inline=True)
    return embed

def gsrs(arg):
    chrome_options = Options()
    chrome_options.add_argument('--headless')  # Run headless Chrome (no UI)
    chrome_options.add_argument('--no-sandbox')
    chrome_options.add_argument('--disable-dev-shm-usage')
    driver = webdriver.Chrome(service=Service(executable_path='/home/spawd/.local/bin/chromedriver'), options=chrome_options)
    try:
        search_url = f'https://gsrs.ncats.nih.gov/ginas/app/beta/browse-substance?search={arg}'
        driver.get(search_url)
        driver.implicitly_wait(10)  # Adjust the wait time as needed
        img_element = driver.find_element(By.CSS_SELECTOR, 'body > app-root > app-base > app-substances-browse > div > div.substance-cards > app-substance-summary-card > mat-card > mat-card-title > a')
        if img_element:
            img_src = img_element.get_attribute('href')
            if img_src:
                stripped = img_src.split('/', -1)[-1:]
                link = f'https://gsrs.ncats.nih.gov/api/v1/substances/render({stripped[0]})?format=png&size=512&stereo=true'
                response = requests.get(link)
                image_bytes = response.content
                image = Image.open(BytesIO(image_bytes))
                image = image.convert('RGBA')
                draw = ImageDraw.Draw(image)
                width, height = image.size
                diagonal = math.sqrt(width**2 + height**2)
                font_size = int(diagonal / 15)
                try:
                    font = ImageFont.truetype('Roboto-Regular.ttf', font_size)
                except IOError:
                    font = ImageFont.load_default()
                bbox = draw.textbbox((0, 0), arg, font=font)
                text_width = bbox[2] - bbox[0]
                text_height = bbox[3] - bbox[1]
                text_x = (width - text_width) / 2
                text_y = (height - text_height) / 2
                watermark_image = Image.new('RGBA', image.size, (0, 0, 0, 0))
                watermark_draw = ImageDraw.Draw(watermark_image)
                watermark_draw.text((text_x, text_y), arg, font=font, fill=(255, 255, 255, 64))
                mask = watermark_image.split()[3]
                image.paste(watermark_image, (0, 0), mask)
                return image
            else:
                return 'No src attribute found in the <img> element'
        else:
            return 'No <img> element found with the specified CSS path'
    finally:
        driver.quit()

def increment_version(config):
    current_version = config['version']
    major, minor, patch = map(int, current_version.split('.'))
    patch += 1
    if patch >= 10:
        patch = 0
        minor += 1
    if minor >= 10:
        minor = 0
        major += 1
    new_version = f'{major}.{minor}.{patch}'
    config['version'] = new_version
    with open(help.path_config_yaml, 'w') as file:
        yaml.dump(config, file)

def invert_colors(image):
    return Image.eval(image, lambda x: 255 - x)

def read_users():
    if not os.path.exists(path_users_yaml):
        raise FileNotFoundError('User file not found.')
    with open(path_users_yaml, 'r') as f:
        return yaml.safe_load(f)

def setup_logging(config: Dict[str, Any]) -> None:
    logging_level = config['logging_level'].upper()
    logging.basicConfig(level=getattr(logging, logging_level))
    if not exists(dirname(path_log)):
        makedirs(dirname(path_log))
    file_handler = logging.FileHandler(path_log)
    file_handler.setLevel(logging_level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger = logging.getLogger()
    logger.setLevel(logging_level)
    logger.addHandler(file_handler)

def unique_pairs(strings_list):
    pairs = list(itertools.combinations(strings_list, 2))
    sorted_pairs = [sorted(list(pair)) for pair in pairs]
    sorted_pairs_overall = sorted(sorted_pairs)
    return sorted_pairs_overall

def validate_generated_code(generated_code):
    forbidden_patterns = [
        'import os', 'import sys', 'class ', 'def setup', 
        'bot.run', 'os.system', 'subprocess', 'eval(', 'exec('
    ]
    for pattern in forbidden_patterns:
        if pattern in generated_code:
            return False
    required_definitions = ['def ']
    if not all(def_ in generated_code for def_ in required_definitions):
        return False
    return True

def write_users(data):
    if not os.path.exists(path_users_yaml):
        os.path.makedirs(os.path.dirname(path_users_yaml))
    with open(path_users_yaml, 'w') as f:
        return yaml.dump(data, f)

