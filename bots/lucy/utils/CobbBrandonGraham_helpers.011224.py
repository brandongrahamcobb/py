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


from bs4 import BeautifulSoup
from collections import defaultdict
from datetime import datetime
from discord import Embed
from discord.ext import commands
from gradio_client import Client
from gtts import gTTS
from io import BytesIO
from matplotlib import pyplot as plt
from openai import AsyncOpenAI
from os import makedirs
from os.path import abspath, dirname, exists, expanduser, isfile, join
from PIL import Image, ImageFont, ImageDraw
from pydub import AudioSegment
from pydub.playback import play
from random import randint
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, DataStructs, Draw, rdDepictor, rdFingerprintGenerator, rdFMCS
rdDepictor.SetPreferCoordGen(True)
from rdkit.Chem.Draw import rdMolDraw2D, SimilarityMaps
from rdkit.DataStructs import FingerprintSimilarity, TanimotoSimilarity
from recipe_scrapers import scrape_me
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from typing import List
from webdriver_manager.chrome import ChromeDriverManager
from typing import List, Optional, Dict, Any
from urllib.parse import urlparse

global logger

import asyncio
import base64
import colorsys
import datetime as dt
import discord
import emoji as emoji_lib
import itertools
import logging
import logging.handlers
import math
import openai
import os
import pubchempy as pcp
import random
import requests
import speech_recognition as sr
import traceback
import unicodedata
import yaml

conversations = defaultdict(list)
current_date = dt.datetime.now().strftime('%d%m%y')

path_home = expanduser('~')
path_config_yaml = join(path_home, '.config', 'vyrtuous', 'config.yaml')
path_log = join(path_home, '.log', 'discord.log')
path_users_yaml = join(path_home, '.config', 'vyrtuous', 'users.yaml')

dir_base = dirname(abspath(__file__))
path_helpers = join(dir_base, 'helpers.py')
path_hybrid = join(dir_base, '..', 'cogs', 'hybrid.py')
path_indica = join(dir_base, '..', 'cogs', 'indica.py')
path_main = join(dir_base, '..', 'main.py')
path_sativa = join(dir_base, '..', 'cogs', 'sativa.py')

helpers_py = load_file(path_helpers)
hybrid_py = load_file(path_hybrid)
indica_py = load_file(path_indica)
main_py = load_file(path_main)
sativa_py = load_file(path_sativa)

def add_watermark(image: BytesIO, watermark_text: str = '~spooky~') -> BytesIO:
    RGB_image = Image.open(image)
    RGBA_image = RGB_image.convert('RGBA')
    draw = ImageDraw.Draw(RGBA_image)
    width, height = RGBA_image.size
    diagonal = math.sqrt(width**2 + height**2)
    font_size = int(diagonal / 15)
    try:
        font = ImageFont.truetype('Roboto-Regular.ttf', font_size)
    except IOError:
        font = ImageFont.load_default()
    while True:
        bbox = draw.textbbox((0, 0), watermark_text, font=font)
        text_width = bbox[2] - bbox[0]
        if text_width <= 512 or font_size <= 1:
            break
        font_size -= 1
        font = ImageFont.truetype('Roboto-Regular.ttf', font_size)
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

def combine(images: list, names: list) -> BytesIO:
    combined_images = []
    for index, (bytes_io, name) in enumerate(zip(images, names)):
        img = Image.open(bytes_io)
        inverted_image = Image.eval(img, lambda x: 255 - x)
        img_bytes = BytesIO()
        inverted_image.save(img_bytes, format='PNG')
        img_bytes_final = add_watermark(img_bytes, watermark_text=name)
        img_final = Image.open(img_bytes_final)
        combined_images.append(img_final)
    widths, heights = zip(*(img.size for img in combined_images))
    total_width = sum(widths)
    max_height = max(heights)
    combined_img = Image.new('RGB', (total_width, max_height))
    x_offset = 0
    for img in combined_images:
        combined_img.paste(img, (x_offset, 0))
        x_offset += img.width
    output = adjust_hue_and_saturation(combined_img, hue_shift=-180, saturation_shift=160)
    output.seek(0)
    return output

def chunk_string(text: str, limit: int = 1800) -> list:
    chunks = []
    while len(text) > limit:
        split_point = text.rfind(' ', 0, limit)
        if split_point == -1:
            split_point = limit
        chunks.append(text[:split_point])
        text = text[split_point:].lstrip()
    if text:
        chunks.append(text)
    return chunks

async def deprecated_create_completion(input_text, sys_input, conversation_id):
    try:
        config = load_config()
        api_key = config['api_keys']['api_key_2']
        ai_client = AsyncOpenAI(api_key=api_key)
        messages = conversations[conversation_id]
        messages.append({'role': 'system', 'content': f'You are Lucy, the Discord bot. Your responses are limited to 2000 characters. Your main.py file is {main_py}. Your cogs are in cogs/ {hybrid_py}, {indica_py}, {sativa_py}. Your helpers are in utils/ {helpers_py}. Your additional System input is: {sys_input}'})
        messages.append({'role': 'user', 'content': input_text})
        stream = await ai_client.chat.completions.create(
            model='gpt-4o-mini',
            messages=messages,
            stream=True
        )
        full_response = ''
        async for chunk in stream:
            content = chunk.choices[0].delta.content
            if content is not None:
                full_response += content
        conversations[conversation_id].append({'role': 'assistant', 'content': full_response})
        yield full_response
    except Exception as e:
        yield traceback.format_exc()

async def create_completion(input_text, conversation_id):
    try:
        config = load_config()
        api_key = config['api_keys']['api_key_2']
        ai_client = AsyncOpenAI(api_key=api_key)
        messages = conversations[conversation_id]
        messages.append({'role': 'user', 'content': input_text})
        stream = await ai_client.chat.completions.create(
            model='o1-mini',
            messages=messages,
            stream=True
        )
        full_response = ''
        async for chunk in stream:
            content = chunk.choices[0].delta.content
            if content is not None:
                full_response += content
        conversations[conversation_id].append({'role': 'assistant', 'content': full_response})
        for i in range(0, len(full_response), 2000):
            yield full_response[i:i + 2000]
    except Exception as e:
        yield traceback.format_exc()

def draw_fingerprint(pair) -> BytesIO:
#    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, countSimulation=True)
#    def get_fp(mol, *args, **kwargs):
#        return mfpgen.GetFingerprint(mol)
    d2d = rdMolDraw2D.MolDraw2DCairo(1024, 1024)
    d2d.prepareMolsBeforeDrawing = False
    Options = d2d.drawOptions()
    Options.prepareMolsBeforeDrawing = False
    Options.includeMetadata = False
    Options.bondLineWidth = 4.0
    d2d.SetDrawOptions(Options)
    mol1 = rdMolDraw2D.PrepareMolForDrawing(pair[0], kekulize=True)
    mol1.UpdatePropertyCache(False)
    mol2 = rdMolDraw2D.PrepareMolForDrawing(pair[1], kekulize=True)
    mol2.UpdatePropertyCache(False)
    fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(mol1, mol2, lambda m, i: SimilarityMaps.GetMorganFingerprint(m, i, radius=2, fpType='bv', nBits=8192), draw2d=d2d, drawingOptions=Options)
#    fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(pair[1], pair[0], get_fp, draw2d=d2d) #colorMap=brighter_color, draw2d=d2d)
    d2d.FinishDrawing()
    drawing = d2d.GetDrawingText()
    output = BytesIO(drawing)
    return output

def draw_watermarked_molecule(molecule) -> BytesIO:
    resolved_name = get_molecule_name(molecule)
    d2d = rdMolDraw2D.MolDraw2DCairo(1024, 1024)
    Options = d2d.drawOptions()
    Options.prepareMolsBeforeDrawing = False
    Options.includeMetadata = False
    Options.bondLineWidth = 4.0
    d2d.SetDrawOptions(Options)
    rdMolDraw2D.SetDarkMode(Options)
    mol = rdMolDraw2D.PrepareMolForDrawing(molecule, kekulize=True)
    mol.UpdatePropertyCache(False)
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    drawing = d2d.GetDrawingText()
    output = add_watermark(BytesIO(drawing), watermark_text=resolved_name)
    return output

def get_proximity(default, input) -> float:
    default_fp = AllChem.GetMorganFingerprintAsBitVect(default, 2)
    input_fp = AllChem.GetMorganFingerprintAsBitVect(input, 2)
    similarity = DataStructs.FingerprintSimilarity(default_fp, input_fp)
    return similarity

def get_mol(arg):
    try:
        compounds = pcp.get_compounds(arg, 'name')
        compound_data = compounds[0].to_dict(properties=['isomeric_smiles'])
        mol = Chem.MolFromSmiles(compound_data['isomeric_smiles'])
        if mol is None:
            raise ValueError('Invalid SMILES string')
        return mol
    except:
        mol = Chem.MolFromSmiles(arg)
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
    return 'Unknown'

def get_recipe(url):
    try:
        scraper = scrape_me(url)
        title = scraper.title()
        total_time = scraper.total_time()
        ingredients = scraper.ingredients()
        instructions = scraper.instructions()
        image = scraper.image()
        servings = scraper.yields()
        source_url = scraper.host()
    except Exception as e:
        return None
    ingredients_formatted = "\n".join([f"• {item}" for item in ingredients])
    instructions_formatted = instructions.replace('\n', '\n')
    embed = Embed(title=title, url=url, color=0x1ABC9C)
    if image:
        embed.set_thumbnail(url=image)
    embed.add_field(name="Servings", value=servings or "N/A", inline=True)
    embed.add_field(name="Total Time", value=f"{total_time} minutes" if total_time else "N/A", inline=True)
    embed.add_field(name="Ingredients", value=ingredients_formatted, inline=False)
    embed.add_field(name="Instructions", value=instructions_formatted, inline=False)
    embed.set_footer(text=f"Source: {source_url}")
    return embed

async def get_other_recipe(url):
    response = requests.get(url)
    if response.status_code != 200:
        return None
    soup = BeautifulSoup(response.content, 'html.parser')
    recipe_div = soup.find('div', id='recipe')
    if not recipe_div:
        return None
    title = recipe_div.find('h2', class_='wprm-recipe-name').text.strip()
    ingredients = []
    ingredient_elements = recipe_div.find_all('ul', class_='wprm-recipe-ingredients')
    for ul in ingredient_elements:
        for li in ul.find_all('li'):
            ingredients.append(li.text.strip())
    instructions = []
    instruction_elements = recipe_div.find_all('ol', class_='wprm-recipe-instructions')
    for ol in instruction_elements:
        for li in ol.find_all('li'):
            instructions.append(li.text.strip())
    embed = discord.Embed(title=title, color=discord.Color.green())
    embed.add_field(name="Ingredients", value='\n'.join(ingredients[:10]), inline=False)
    embed.add_field(name="Instructions", value='\n'.join(instructions[:5]), inline=False)
    return embed


def get_scripture(version: str, reference: str):
     try:
         BIBLE_IDS = {
             'esv': 'de4e12af7f28f599-02',
             'nkjv': 'de4e12af7f28f599-01',
             'niv': '06125adad2d5898a-01',
         }
         version = version.lower()
         if version in BIBLE_IDS:
             bible_id = BIBLE_IDS[version]
             api = f'https://api.scripture.api.bible/v1/bibles/{bible_id}/search?query={reference}'
             response = requests.get(api, headers=get_scripture_headers())
             if response.ok:
                 json = response.json()
                 passages = json.get('data', {}).get('passages', [])
                 soup = BeautifulSoup(passages[0].get('content'), 'html.parser')
                 soup.get_text()
                 cleaned_content = soup.get_text()
                 message = f'**{reference}** ({version.upper()})\n{cleaned_content}'
                 return message
         else:
             response = requests.get(f'https://api.alquran.cloud/v1/ayah/{reference}/en.asad', headers=get_scripture_headers())
             if response.ok:
                 json = response.json()
                 message = f"**{reference}** ({version.upper()})\n{json['data']['text']}"
                 return message
     except Exception as e:
         return e

def get_scripture_headers():
        headers = {
            'User-Agent': 'spooky',
            'api-key': '2eb327f99245cd3d68da55370656d6e2'
        }
        return headers

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

async def join_voice_channel(channel_id):
    channel = self.bot.get_channel(channel_id)
    if isinstance(channel, discord.VoiceChannel):
        return await channel.connect()
    else:
        raise Exception("Channel is not a voice channel.")

def google(query: str, num_results: int = 5):
    headers = {
        "User-Agent": (
            "Chrome/91.0.4472.124 Arch Linux 6.12.1-arch1-1"
        )
    }
    search_url = "https://www.google.com/search"
    params = {"q": query, "num": num_results}

    try:
        response = requests.get(search_url, headers=headers, params=params)
        response.raise_for_status()  # Raise an error for bad HTTP responses
        soup = BeautifulSoup(response.text, "html.parser")

        results = []
        for g in soup.find_all("div", class_="tF2Cxc"):
            title = g.find("h3").text if g.find("h3") else "No title"
            link = g.find("a")["href"] if g.find("a") else "No link"
            results.append({"title": title, "link": link})
            if len(results) >= num_results:
                break

        return results
    except requests.exceptions.RequestException as e:
        print(f"Error during the web request: {e}")
        return []

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

def load_config():
    if exists (path_config_yaml):
        with open(path_config_yaml, 'r') as f:
            data = yaml.safe_load(f)
            if data:
                return data
    else:
        return None

def load_file(path_to_file):
    if not exists(path_to_file):
        raise FileNotFoundError(f"The file at '{path_to_file}' does not exist.")
    try:
        with open(path_to_file, 'r', encoding='utf-8') as file:
            content = file.read()
        return content
    except Exception as e:
        raise IOError(f"An error occurred while reading the file: {e}")

def load_users():
    if exists(path_users_yaml):
        with open(join(path_home, 'Downloads', 'users.yaml'), 'r') as f:
            data = yaml.safe_load(f)
            if data:
                return set(data.get('users', []))
    else:
        return None

def increment_version(config: Dict[str, Any]):
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
    with open(path_config_yaml, 'w') as file:
        yaml.dump(config, file)

def save_users(data):
    with open(join(path_home, 'Downloads', 'users.yaml'), 'w') as f:
        yaml.dump({'users': list(data)}, f)

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

def stable_cascade(prompt):
    try:
        client = Client('multimodalart/stable-cascade')
        result = client.predict(
            prompt=prompt,
            negative_prompt='',
            seed=randint(0, 2147483647),
            width=1024,
            height=1024,
            prior_num_inference_steps=20,
            prior_guidance_scale=4,
            decoder_num_inference_steps=10,
            decoder_guidance_scale=0,
            num_images_per_prompt=1,
            api_name="/run",
        )
        return discord.File(result, 'image.webp')
    except ConnectionError as conn_err:
        print(f"Connection error: {conn_err}")
        return "Failed to connect to the server. Please try again later."
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return f"An error occurred: {e}"

async def text_to_voice(self, text, voice_client):
    audio_file = 'response.mp3'
    tts = gTTS(text=text, lang='en', slow=False)
    tts.save(audio_file)
    if voice_client.is_connected():
        audio_source = discord.FFmpegPCMAudio(audio_file)
        voice_client.play(audio_source, after=lambda e: os.remove(audio_file) if e is None else None)
    else:
        await voice_client.disconnect()

def unique_pairs(strings_list):
    pairs = list(itertools.combinations(strings_list, 2))
    sorted_pairs = [sorted(list(pair)) for pair in pairs]
    sorted_pairs_overall = sorted(sorted_pairs)
    return sorted_pairs_overall

async def voice_to_text(ctx):
    r = sr.Recognizer()
    with sr.Microphone() as source:
        r.adjust_for_ambient_noise(source)
        print("Listening...")
        try:
            audio = r.listen(source, timeout=30)
            text = r.recognize_google(audio)
            await ctx.send(f"You said: {text}")
            return text
        except sr.UnknownValueError:
            await ctx.send("I could not understand the audio input.")
            return None
        except Exception as e:
            print(f"Error: {e}")
            await ctx.send("An error occurred while processing the audio.")
            return None
