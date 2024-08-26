""" user_cog.py
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

from bot.utils.helpers import get_emoji
from bot.utils.helpers import get_sds
from bot.utils.helpers import search
from discord.ext import commands
from googletrans import Translator, LANGUAGES
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from webdriver_manager.chrome import ChromeDriverManager
from PIL import Image, ImageFont, ImageDraw

import discord
import json
import io
import math
import PIL
import pubchempy as pcp
import requests
import traceback

class UserCog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.translator = Translator()
        self.user_translation_preferences = {}

    @commands.Cog.listener()
    async def on_message(self, message: discord.Message):
        if message.author.bot:
            return
        prefs = self.user_translation_preferences.get(message.author.id)
        if prefs:
            target_lang_code, source_lang_code = prefs
            try:
                translated = self.translator.translate(message.content, dest=target_lang_code, src=source_lang_code)
                await message.channel.send(f'{translated.text}')
            except Exception as e:
                await message.channel.send(f'Error translating message: {e}')

    def get_language_code(self, language_name):
        language_name = language_name.lower()
        for lang_code, lang_name in LANGUAGES.items():
            if lang_name.lower() == language_name:
                return lang_code
        return None

    def gsrs(self, query):
        chrome_options = Options()
        chrome_options.add_argument("--headless")  # Run headless Chrome (no UI)
        chrome_options.add_argument("--no-sandbox")
        chrome_options.add_argument("--disable-dev-shm-usage")
        driver = webdriver.Chrome(service=Service(executable_path="/home/spawd/.local/bin/chromedriver"), options=chrome_options)
        try:
            search_url = f"https://gsrs.ncats.nih.gov/ginas/app/beta/browse-substance?search={query}"
            driver.get(search_url)
            driver.implicitly_wait(10)  # Adjust the wait time as needed
            img_element = driver.find_element(By.CSS_SELECTOR, "body > app-root > app-base > app-substances-browse > div > div.substance-cards > app-substance-summary-card > mat-card > mat-card-title > a")
            if img_element:
                img_src = img_element.get_attribute('href')
                if img_src:
                    stripped = img_src.split('/', -1)[-1:]
                    link = f'https://gsrs.ncats.nih.gov/api/v1/substances/render({stripped[0]})?format=png&size=512&stereo=true'
                    response = requests.get(link)
                    image_bytes = response.content
                    image = Image.open(io.BytesIO(image_bytes))
                    image = image.convert("RGBA")
                    draw = ImageDraw.Draw(image)
                    width, height = image.size
                    diagonal = math.sqrt(width**2 + height**2)
                    font_size = int(diagonal / 15)  # Adjust this factor as needed
                    try:
                        font = ImageFont.truetype("arial.ttf", font_size)
                    except IOError:
                        font = ImageFont.load_default(font_size)
                    bbox = draw.textbbox((0, 0), query, font=font)
                    text_width = bbox[2] - bbox[0]
                    text_height = bbox[3] - bbox[1]
                    text_x = (width - text_width) / 2
                    text_y = (height - text_height) / 2
                    watermark_image = Image.new("RGBA", image.size, (0, 0, 0, 0))
                    watermark_draw = ImageDraw.Draw(watermark_image)
                    watermark_draw.text((text_x, text_y), query, font=font, fill=(255, 255, 255, 64))
                    mask = watermark_image.split()[3]
                    image.paste(watermark_image, (0, 0), mask)
                    return image
                else:
                    return "No src attribute found in the <img> element"
            else:
                return "No <img> element found with the specified CSS path"
        finally:
            driver.quit()

    @commands.command(name='get')
    async def get(self, ctx, *, query: str):
        results = search(query, 'web')
        embed = discord.Embed()
        for result in results:
            embed.add_field(name=result['title'], value=result['link'], inline=False)
        await ctx.send(embed=embed)

    @commands.command(name='info')
    async def info(self, ctx: commands.Context, argument: str = None):
        if argument is None:
            await ctx.send('The command is used !info <MOLECULE> or !info :emoji:')
            return
        if isinstance(argument, discord.Emoji):
            try:
                await ctx.send(get_emoji(argument))
            except Exception as e:
                await ctx.send(e)
        else:
            try:
                await ctx.send(get_sds(argument))
            except Exception as e:
                await ctx.send(e)

    @commands.command()
    async def languages(self, ctx):
        supported_languages = ', '.join(LANGUAGES.values())
        await ctx.send(f'Supported languages are:\n{supported_languages}')

    @commands.command(name='search')
    async def search(self, ctx, *, query=None):
        try:
            try:
                watermarked_image = self.gsrs(query)
                with io.BytesIO() as image_binary:
                    watermarked_image.save(image_binary, format='PNG')
                    image_binary.seek(0)
                    await ctx.send(file=discord.File(fp=image_binary, filename='watermarked_image.png'))
            except:
                await ctx.send(traceback.format_exc())
        except:
            if query is None:
                return
            await ctx.send(f'{query} is an unknown molecule.')

    @commands.command(name='smiles')
    async def smiles(self, ctx: commands.Context, *, chemical_name: str):
        compounds = pcp.get_compounds(chemical_name, 'name')
        compound = compounds[0]
        isomeric_smiles = compound.isomeric_smiles
        await ctx.send(f'The isomeric SMILES for {chemical_name} is: {isomeric_smiles}')

    @commands.command()
    async def translate(self, ctx, toggle: str, target_lang: str = 'english', source_lang: str = 'auto'):
        if toggle.lower() == 'on':
            target_lang_code = self.get_language_code(target_lang)
            source_lang_code = self.get_language_code(source_lang)
            if target_lang_code is None or source_lang_code is None:
                await ctx.send(f'{ctx.author.mention}, please specify valid language names.')
                return
            self.user_translation_preferences[ctx.author.id] = (target_lang_code, source_lang_code)
            await ctx.send(f'{ctx.author.mention}, translation enabled from {source_lang} to {target_lang}.')
        elif toggle.lower() == 'off':
            self.user_translation_preferences[ctx.author.id] = None
            await ctx.send(f'{ctx.author.mention}, translation disabled.')
        else:
            await ctx.send(f'{ctx.author.mention}, please specify "on" or "off".')

async def setup(bot: commands.Bot):
    await bot.add_cog(UserCog(bot))
