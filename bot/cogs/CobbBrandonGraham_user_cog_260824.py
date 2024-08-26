''' user_cog.py
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
from bot.utils.helpers import get_images
from bot.utils.helpers import get_emoji
from bot.utils.helpers import get_mol
from bot.utils.helpers import get_sds
from bot.utils.helpers import read_config
from bot.utils.helpers import read_users
from bot.utils.helpers import write_users
from bot.utils.helpers import off_peak_hours_required
from bot.utils.helpers import unique_pairs
from discord.ext import commands
from functools import wraps
from googleapiclient.discovery import build
from googletrans import Translator, LANGUAGES
from io import BytesIO
from PIL import Image, ImageFont, ImageDraw
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, DataStructs, Draw
from rdkit.Chem.Draw import rdMolDraw2D, SimilarityMaps
from rdkit.DataStructs import FingerprintSimilarity, TanimotoSimilarity
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import IPythonConsole
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from typing import List
from webdriver_manager.chrome import ChromeDriverManager

import colorsys
import discord
import io
import json
import math
import os
import PIL
import pubchempy as pcp
import random
import rdkit
import requests
import tempfile
import traceback
import yaml

class UserCog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.channels = bot.get_all_channels()
        self.country_references = {
            'USA': get_mol('LSD'),
        }
        self.guild = self.bot.get_guild(self.bot.testing_guild_id)
#        self.headers = read_config()['headers']
        self.tags = {}
        self.translator = Translator()
        self.user_translation_preferences = {}
        self.user_id = '154749533429956608'
#        self.user_data = read_users()
        self.stacks = {}
        self.lysergic_acid_diethylamide = get_mol('LSD')

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

    def add_watermark(self, image: BytesIO, watermark_text: str = '~LUCY~') -> BytesIO:
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

    def adjust_hue_and_saturation(self, image, hue_shift, saturation_shift) -> BytesIO:
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

    def combine(self, bytes1: BytesIO, bytes2: BytesIO) -> BytesIO:
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
        inverted_image = self.invert_colors(combined_img)
        output = self.adjust_hue_and_saturation(inverted_image, hue_shift=-180, saturation_shift=160)
        output.seek(0)
        return output

    def draw_fingerprint(self, pair) -> BytesIO:
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

    def draw_watermarked_molecule(self, molecule) -> BytesIO:
        resolved_name = self.get_molecule_name(molecule)
        d2d = rdMolDraw2D.MolDraw2DCairo(1024, 1024)
        Options = d2d.drawOptions()
        d2d.SetDrawOptions(Options)
        rdMolDraw2D.SetDarkMode(Options)
        mol = rdMolDraw2D.PrepareMolForDrawing(molecule, kekulize=True)
        mol.UpdatePropertyCache(False)
        d2d.DrawMolecule(mol)
        d2d.FinishDrawing()
        drawing = d2d.GetDrawingText()
        output = self.add_watermark(BytesIO(drawing), watermark_text=resolved_name)
        return output

    def get_language_code(self, language_name):
        language_name = language_name.lower()
        for lang_code, lang_name in LANGUAGES.items():
            if lang_name.lower() == language_name:
                return lang_code
        return None

    def get_proximity(self, default, input) -> float:
        default_fp = AllChem.GetMorganFingerprintAsBitVect(default, 2)
        input_fp = AllChem.GetMorganFingerprintAsBitVect(input, 2)
        similarity = DataStructs.FingerprintSimilarity(default_fp, input_fp)
        return similarity

    def get_molecule_name(self, molecule) -> str:
        smiles = Chem.MolToSmiles(molecule)
        if smiles:
            compounds = pcp.get_compounds(smiles, 'smiles') #, record_type='3d')
            compound_data = compounds[0].to_dict(properties=['synonyms'])
            if not compounds:
                raise ValueError('No compound found for the given SMILES string')
            return compound_data['synonyms'][0]

    def gsrs(self, arg):
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
                    image = Image.open(io.BytesIO(image_bytes))
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

    def invert_colors(self, image):
        return Image.eval(image, lambda x: 255 - x)

    @commands.command(name='covid', help='Get the latest COVID-19 statistics globally or for a specific country')
    async def covid(self, ctx, country: str = 'global'):
        '''Fetches COVID-19 statistics from the disease.sh API.'''
        if country.lower() == 'global':
            url = 'https://disease.sh/v3/covid-19/all'
        else:
            url = f'https://disease.sh/v3/covid-19/countries/{country}'

        response = requests.get(url)
        data = response.json()

        if response.status_code == 200:
            embed = discord.Embed(title=f'COVID-19 Stats for {country.title()}', color=discord.Color.blue())
            embed.add_field(name='Cases', value=f'{data['cases']:,}', inline=True)
            embed.add_field(name='Deaths', value=f'{data['deaths']:,}', inline=True)
            embed.add_field(name='Recovered', value=f'{data['recovered']:,}', inline=True)
            embed.add_field(name='Active', value=f'{data['active']:,}', inline=True)
            embed.add_field(name='Critical', value=f'{data['critical']:,}', inline=True)
            embed.add_field(name='Today\'s Cases', value=f'{data['todayCases']:,}', inline=True)
            embed.add_field(name='Today\'s Deaths', value=f'{data['todayDeaths']:,}', inline=True)
            await ctx.send(embed=embed)
        else:
            await ctx.send('Failed to retrieve data. Please try again later or check the country name.')

#    @commands.cooldown(rate=3, per=1.0, type=commands.BucketType.channel)
    @commands.command(name='draw')
    @off_peak_hours_required()
    async def draw(self, ctx, *args) -> None:
 #       if await self.on_cooldown(ctx, 'draw'):
  #          return
        async with ctx.typing():
            if len(args) == 1:
                 image = self.draw_watermarked_molecule(get_mol(args[0]))
                 await ctx.send(file=discord.File(image, f'{args[0]}.png'))
            pairs = unique_pairs(args)
            for pair in pairs:
                 mol = get_mol(pair[0])
                 refmol = get_mol(pair[1])
                 image = self.draw_watermarked_molecule(mol)
                 proximity = self.get_proximity(default=mol, input=refmol)
                 await ctx.send(f'The Tanimoto Similarity of {self.get_molecule_name(mol)} and {self.get_molecule_name(refmol)} is {proximity:.2f}.')
                 fingerprints = []
                 fingerprints.append(
                     self.draw_fingerprint(
                         [mol, refmol]
                     )
                 )
                 fingerprints.append(
                     self.draw_fingerprint(
                         [refmol, mol]
                     )
                 )
                 combined_image = self.combine(fingerprints[0], fingerprints[1])
                 await ctx.send(file=discord.File(combined_image, f'molecule_comparison.png'))

    @commands.command(name='get', hidden=True)
    @commands.is_owner()
    async def get(self, ctx, *, arg: str):
        result = get_images(arg)
        await ctx.send(result)

    @commands.command(name='info', hidden=True)
    async def info(self, ctx: commands.Context, argument):
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
                await ctx.send(embed=get_sds(argument))
            except Exception as e:
                await ctx.send(e)

    @commands.command(name='languages', hidden=True)
    @commands.is_owner()
    async def languages(self, ctx):
        supported_languages = ', '.join(LANGUAGES.values())
        await ctx.send(f'Supported languages are:\n{supported_languages}')

    @commands.command(name='logp', hidden=True)
    @commands.is_owner()
    async def logp(self, ctx, *args):
        try:
            for arg in args:
                compounds = pcp.get_compounds(arg, 'name')
                compound_data = compounds[0].to_dict(properties=['isomeric_smiles'])
                mol = Chem.MolFromSmiles(compound_data['isomeric_smiles'])
                log_p = Crippen.MolLogP(mol)
                await ctx.send(f'Your octanol:water coefficient is: {log_p}')
        except Exception as e:
            await ctx.send(f'Error fetching data: {e}')

    @commands.command(name='search', hidden=True)
    async def search(self, ctx, *, arg=None):
        try:
            try:
                watermarked_image = self.gsrs(arg)
                with io.BytesIO() as image_binary:
                    watermarked_image.save(image_binary, format='PNG')
                    image_binary.seek(0)
                    await ctx.send(file=discord.File(fp=image_binary, filename='watermarked_image.png'))
            except:
                await ctx.send(traceback.format_exc())
        except:
            if arg is None:
                return
            await ctx.send(f'{arg} is an unknown molecule.')

    @commands.command(name='smiles', hidden=True)
    async def smiles(self, ctx: commands.Context, *, chemical_name: str):
        compounds = pcp.get_compounds(chemical_name, 'name')
        compound = compounds[0]
        isomeric_smiles = compound.isomeric_smiles
        await ctx.send(f'The isomeric SMILES for {chemical_name} is: {isomeric_smiles}')

    def get_user_stack(self, user_id: int):
        return self.stacks.get(user_id, [])

    def set_user_stack(self, user_id: int, molecule_names: List[str]):
        mol_objects = [get_mol(name) for name in molecule_names]
        self.stacks[user_id] = mol_objects

    @commands.command(name='stack', description='Submit one molecule to test or multiple to set a stack.')
    async def stack(self, ctx: commands.Context, *molecules):
        if not molecules:
            await ctx.send(f'Current negative node: LSD')
            await ctx.send([self.get_molecule_name(mol) for mol in self.get_user_stack(ctx.author.id)])
        if ctx.author.id in self.stacks and len(molecules) == 1:
            eximity = self.get_proximity(get_mol(molecules[0]), self.lysergic_acid_diethylamide)
            defender = max(self.get_user_stack(ctx.author.id), key=lambda mol: self.get_proximity(mol, get_mol(molecules[0])))
            proximity = self.get_proximity(get_mol(molecules[0]), defender)
            if eximity > proximity:
                await ctx.send(
                    f'{self.get_molecule_name(defender)} is closest ({proximity:.3f} vs. {eximity:.3f}) but, {molecules[0]} is \N{REVERSED HAND WITH MIDDLE FINGER EXTENDED} for your stack.'
                )
            else:
                await ctx.send(
                    f'Before you take {molecules[0]}, contemplate how it interacts with {self.get_molecule_name(defender)}! \N{OK HAND SIGN}\n'
                    f'Please review resources on this compound with !info {molecules[0]}'
                )
        elif len(molecules) > 1:
            self.set_user_stack(user_id=ctx.author.id, molecule_names=molecules)
            await ctx.send('Stack overwritten.')

    @commands.command(name='translate', hidden=True)
#    @commands.is_owner()
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
            await ctx.send(f'{ctx.author.mention}, please specify `on` or `off`.')


async def setup(bot: commands.Bot):
    await bot.add_cog(UserCog(bot))
