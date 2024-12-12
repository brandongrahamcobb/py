''' helpers.py  The purpose of this program is to provide generic python methods.
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
from datetime import datetime
from discord import Embed
from gradio_client import Client
from io import BytesIO
from os import makedirs
from os.path import abspath, dirname, exists, expanduser, isfile, join
from PIL import Image, ImageFont, ImageDraw
from recipe_scrapers import scrape_me
from typing import Dict, Any
from urllib.parse import urlparse

global logger

import colorsys
import datetime as dt
import discord
import emoji as emoji_lib
import itertools
import json
import logging
import logging.handlers
import math
import os
from random import randint
import requests
import traceback
import unicodedata
import yaml

current_date = dt.datetime.now().strftime('%d%m%y')

path_home = expanduser('~')
path_config_yaml = join(path_home, '.config', 'vyrtuous', 'config.yaml')
path_log = join(path_home, '.log', 'discord.log')

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

def append_to_jsonl(data, filename: str) -> None:
    """Append a json payload to the end of a jsonl file."""
    json_string = json.dumps(data)
    with open(filename, "a") as f:
        f.write(json_string + "\n")

def chunk_string(text: str, limit: int = 1800) -> list:
    chunks = []
    while len(text) > limit:
        split_point = text.rfind('\n', 0, limit)
        if split_point == -1:
            split_point = limit
        chunks.append(text[:split_point])
        text = text[split_point:].lstrip()
    if text:
        chunks.append(text)
    return chunks

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

def create_embed(title, description, color=0x00ff00):
    from discord import Embed
    return Embed(title=title, description=description, color=color)

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

def google(query: str, num_results: int = 5):
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
    }
    search_url = "https://www.google.com/search"
    params = {"q": query, "num": num_results}
    try:
        response = requests.get(search_url, headers=headers, params=params)
        response.raise_for_status()  # Raise an error for bad HTTP responses
        soup = BeautifulSoup(response.text, "html.parser")
        results = []
        for g in soup.find_all("div", class_="g"):  # Updated class name
            title = g.find("h3").text if g.find("h3") else "No title"
            link = g.find("a")["href"] if g.find("a") else "No link"
            results.append({"title": title, "link": link})
            if len(results) >= num_results:
                break
        return results
    except requests.exceptions.RequestException as e:
        print(f"Error during the web request: {e}")
        return []

def load_contents(path_to_file):
    if not exists(path_to_file):
        raise FileNotFoundError(f"The file at '{path_to_file}' does not exist.")
    try:
        with open(path_to_file, 'r', encoding='utf-8') as file:
            content = file.read()
        return content
    except Exception as e:
        raise IOError(f"An error occurred while reading the file: {e}")

def load_json(path_to_file):
    if not os.path.exists(path_to_file):
        return {}
    with open(path_to_file, 'r') as f:
        return json.load(f) or {}

def load_yaml(path_to_file):
    if not os.path.exists(path_to_file):
        return {}
    with open(path_to_file, 'r') as f:
        return yaml.safe_load(f) or {}

def increment_infraction(user_id, users):
    users[user_id] = users.get(user_id, 0) + 1
    return users[user_id]

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

async def save_json(path_to_file, data):
    with open(path_to_file, 'w') as f:
        json.dump(data, f)

async def save_yaml(path_to_file, data):
    with open(path_to_file, 'w') as f:
        yaml.dump(data, f)

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

def unique_pairs(strings_list):
    pairs = list(itertools.combinations(strings_list, 2))
    sorted_pairs = [sorted(list(pair)) for pair in pairs]
    sorted_pairs_overall = sorted(sorted_pairs)
    return sorted_pairs_overall
