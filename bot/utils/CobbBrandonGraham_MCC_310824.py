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

from bot.main import CustomBot
from datetime import datetime
from discord.ext import commands
from googleapiclient.discovery import build
from io import BytesIO
from os import makedirs
from os.path import abspath, dirname, exists, expanduser, isfile, join
from PIL import Image, ImageFont, ImageDraw
from typing import List
from typing import List, Optional, Dict, Any

global logger

import asyncio
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
import unicodedata
import yaml

config = CustomBot._get_config()

current_date = dt.datetime.now().strftime('%d%m%y')

home = expanduser('~')

path_base = join(home, 'Documents', 'src', 'mcc')
path_ai_cog = join(path_base, 'cogs', 'ai_cog.py')
path_config_yaml = join(home, '.config', 'MCC', 'config.yaml')
path_log = join(home, '.log', 'discord.log')
path_users_yaml = join(home, '.config', 'MCC', 'users.yaml')

def add_watermark(image: BytesIO, watermark_text: str = '~MCC STEMgineering Club~') -> BytesIO:
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

async def handle_command_error(ctx: commands.Context, error: commands.CommandError):
    log_id = str(uuid.uuid4())
    logger.error(f"Log ID {log_id} - Error in command {ctx.command}: {error}")
    logger.error(f"Full stack trace:\n{traceback.format_exc()}")
    await ctx.send(f"An unexpected error occurred (Log ID: `{log_id}`). Please report this to the support team.")

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

def write_users(data: Dict[str, Any]):
    if not os.path.exists(path_users_yaml):
        os.path.makedirs(os.path.dirname(path_users_yaml))
    with open(path_users_yaml, 'w') as f:
        return yaml.dump(data, f)

