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

from bot.utils.helpers import calculate_max_font_size
from bot.utils.helpers import get_emoji
from bot.utils.helpers import get_sds
from bot.utils.helpers import get_version
from bot.utils.helpers import load_config
from bot.utils.helpers import search
from discord.ext import commands
from googletrans import Translator, LANGUAGES
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import SimilarityMaps
from rdkit.DataStructs import FingerprintSimilarity, TanimotoSimilarity
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from webdriver_manager.chrome import ChromeDriverManager

import asyncio
import discord
import json
import io
import pubchempy as pcp
import math
import random
import requests

class UserCog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.translator = Translator()
        self.user_translation_preferences = {}

    def monograph(molecule: str) -> io.BytesIO():
        d2d = Draw.MolDraw2DCairo(512, 512)
        compounds = pcp.get_compounds(molecule, 'name')
        compound_data = compounds[0].to_dict(properties=['isomeric_smiles'])
        text = compounds[0].synonyms[0]
        mol = Chem.MolFromSmiles(compound_data['isomeric_smiles'])
        fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(mol, mol, lambda m, i: SimilarityMaps.GetMorganFingerprint(m, i, radius=2, fpType='bv', nBits=8192), draw2d=d2d)
        d2d.FinishDrawing()
        buf = io.BytesIO(d2d.GetDrawingText())
        image1 = Image.open(buf)
        image1 = image1.convert('RGBA')
        width, height = image1.size
        draw = ImageDraw.Draw(image1)
#        for x in range(image1.width):
 #           for y in range(image1.height):
  #              pixel = image1.getpixel((x, y))
   #             if sum(pixel) >= 1020:
    #                draw.point((x, y), fill=(255, 255, 255, 128))
        font_size = min(width, height) // 10  # starting font size (adjust as needed)
        font = ImageFont.load_default(font_size)
        text_bbox = draw.textbbox((0, 0), text, font=font)
        text_width = text_bbox[2] - text_bbox[0]
        text_height = text_bbox[3] - text_bbox[1]
        while text_width > width - 20:
            font_size -= 1
            font = ImageFont.load_default(font_size)
            text_bbox = draw.textbbox((0, 0), text, font=font)
            text_width = text_bbox[2] - text_bbox[0]
            text_height = text_bbox[3] - text_bbox[1]
        text_x = ((width - text_width) / 2)
        text_y = height - text_height - 60
        draw.text((text_x, text_y), text, font=font, fill=(152,251,203))
        draw.text((text_x + 2, text_y - 2), text, font=font, fill=(152,251,203))
        draw.text((text_x + 4, text_y - 4), text, font=font, fill=(152,251,203))
        output_buffer = io.BytesIO()
        image1.save(output_buffer, format='PNG')
        output_buffer.seek(0)
        return output_buffer

    def to_smiles(molecule: str):
        mol = Chem.MolFromSmiles(molecule)
        img = Draw.MolToImage(mol, size=(512, 512))
        draw = ImageDraw.Draw(img)
        font = ImageFont.load_default(40)
        label_1 = 'Unknown'
        width, height = img.size
        for x in range(img.width):
            for y in range(img.height):
                pixel = img.getpixel((x, y))
                if sum(pixel) >= 1020:
                    draw.point((x, y), fill=(0, 0, 0, 0))
        draw.text(((width - draw.textlength(label_1, font=font)) / 2, height - 60), label_1, fill=(152,251,203), font=font)
        img_bytes = io.BytesIO()
        img.save(img_bytes, format='PNG')
        img_bytes.seek(0)
        return img_bytes

    def monograph_three(args):
        try:
            d2d1 = Draw.MolDraw2DCairo(512, 512)
            d2d2 = Draw.MolDraw2DCairo(512, 512)
            compounds = pcp.get_compounds(args[0], 'name')
            if compounds:
                label_1 = compounds[0].synonyms[0]
                compound_data = compounds[0].to_dict(properties=['isomeric_smiles'])
                mol = Chem.MolFromSmiles(compound_data['isomeric_smiles'])
            else:
                label_1 = 'Unknown'
                mol = Chem.MolFromSmiles(args[0])
            ref_compounds = pcp.get_compounds(args[1], 'name')
            if ref_compounds:
                label_2 = ref_compounds[0].synonyms[0]
                ref_compound_data = ref_compounds[0].to_dict(properties=['isomeric_smiles'])
                refmol = Chem.MolFromSmiles(ref_compound_data['isomeric_smiles'])
            else:
                label_2 = 'Unknown'
                refmol = Chem.MolFromSmiles(args[1])
            fp1 = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            fp2 = AllChem.GetMorganFingerprintAsBitVect(refmol, 2, nBits=2048)
            similarity = TanimotoSimilarity(fp1, fp2)
            fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(refmol, mol, lambda m, i: SimilarityMaps.GetMorganFingerprint(m, i, radius=2, fpType='bv', nBits=8192), draw2d=d2d1)
            d2d1.FinishDrawing()
            buf = io.BytesIO(d2d1.GetDrawingText())
            image1 = Image.open(buf)
            image1 = image1.convert('RGBA')
            fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(mol, refmol, lambda m, i: SimilarityMaps.GetMorganFingerprint(m, i, radius=2, fpType='bv', nBits=8192), draw2d=d2d2)
            d2d2.FinishDrawing()
            buf = io.BytesIO(d2d2.GetDrawingText())
            image2 = Image.open(buf)
            image2 = image2.convert('RGBA')
            width, height = image2.size
            new_image = Image.new('RGBA', (width + width, height), (0, 0, 0, 0))
            new_image.paste(image1, (0, 0), image1)
            new_image.paste(image2, (width, 0), image2)
            original = new_image
            transparent_version = original.copy()
            alpha = transparent_version.split()[3]
            alpha = alpha.point(lambda p: min(p, 128))
            transparent_version.putalpha(alpha)
            combined = Image.new('RGBA', original.size)
            combined.paste(transparent_version, (10, 10), mask=transparent_version)
            combined.paste(original, (0, 0))
            draw = ImageDraw.Draw(combined)
            font_size = calculate_max_font_size(label_1, 512)
            font = ImageFont.load_default(font_size)
            draw.text(((width - draw.textlength(label_1, font=font)) / 2, 10), label_1, fill='gray', font=font)
            font_size = calculate_max_font_size(label_2, 512)
            font = ImageFont.load_default(font_size)
            draw.text((width + (width - draw.textlength(label_2, font=font)) / 2, 10), label_2, fill='gray', font=font)
            output_buffer = UserCog.photoshop(draw, combined, f'{(similarity * 100):.0f}%')
            return output_buffer, similarity
        except Exception as e:
            print(e)

    def monograph_two(args):
        try:
            d2d1 = Draw.MolDraw2DCairo(512, 512)
            d2d2 = Draw.MolDraw2DCairo(512, 512)
            compounds = pcp.get_compounds(args[0], 'name')
            if compounds:
                label_1 = compounds[0].synonyms[0]
                compound_data = compounds[0].to_dict(properties=['isomeric_smiles'])
                mol = Chem.MolFromSmiles(compound_data['isomeric_smiles'])
            else:
                label_1 = 'Unknown'
                mol = Chem.MolFromSmiles(args[0])
            ref_compounds = pcp.get_compounds(args[1], 'name')
            if ref_compounds:
                label_2 = ref_compounds[0].synonyms[0]
                ref_compound_data = ref_compounds[0].to_dict(properties=['isomeric_smiles'])
                refmol = Chem.MolFromSmiles(ref_compound_data['isomeric_smiles'])
            else:
                label_2 = 'Unknown'
                refmol = Chem.MolFromSmiles(args[1])
            fp1 = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            fp2 = AllChem.GetMorganFingerprintAsBitVect(refmol, 2, nBits=2048)
            similarity = TanimotoSimilarity(fp1, fp2)
            fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(refmol, mol, lambda m, i: SimilarityMaps.GetMorganFingerprint(m, i, radius=2, fpType='bv', nBits=8192), draw2d=d2d1)
            d2d1.FinishDrawing()
            buf = io.BytesIO(d2d1.GetDrawingText())
            image1 = Image.open(buf)
            image1 = image1.convert('RGBA')
            fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(mol, refmol, lambda m, i: SimilarityMaps.GetMorganFingerprint(m, i, radius=2, fpType='bv', nBits=8192), draw2d=d2d2)
            d2d2.FinishDrawing()
            buf = io.BytesIO(d2d2.GetDrawingText())
            image2 = Image.open(buf)
            image2 = image2.convert('RGBA')
            width, height = image2.size
            new_image = Image.new('RGBA', (width + width, height), (0, 0, 0, 0))
            new_image.paste(image1, (0, 0), image1)
            new_image.paste(image2, (width, 0), image2)
            original = new_image
            transparent_version = original.copy()
            alpha = transparent_version.split()[3]
            alpha = alpha.point(lambda p: min(p, 128))
            transparent_version.putalpha(alpha)
            combined = Image.new('RGBA', original.size)
            combined.paste(transparent_version, (10, 10), mask=transparent_version)
            combined.paste(original, (0, 0))
            draw = ImageDraw.Draw(combined)
            font_size = calculate_max_font_size(label_1, 512)
            font = ImageFont.load_default(font_size)
            draw.text(((width - draw.textlength(label_1, font=font)) / 2, combined.height - 60), label_1, fill='gray', font=font)
            font_size = calculate_max_font_size(label_2, 512)
            font = ImageFont.load_default(font_size)
            draw.text((width + (width - draw.textlength(label_2, font=font)) / 2, combined.height - 60), label_2, fill='gray', font=font)
            output_buffer = UserCog.photoshop(draw, combined, f'{(similarity * 100):.0f}%')
            return output_buffer, similarity
        except Exception as e:
            print(e)

    def photoshop(draw, image, text):
        font = ImageFont.load_default(40)
        text_bbox = draw.textbbox((0, 0), text, font=font)
        text_width, text_height = text_bbox[2] - text_bbox[0], text_bbox[3] - text_bbox[1]
        text_x = (image.width - text_width) // 2
        text_y = (image.height - text_height) // 2
        draw.text((text_x, text_y), text, fill='gray', font=font)
        output_buffer = io.BytesIO()
        image.save(output_buffer, format='PNG')
        output_buffer.seek(0)
        return output_buffer

    @commands.command(description='Compare molecules, fetch molecules.')
    async def draw(self, ctx: commands.Context, *args) -> None:
        if len(args) == 1 or len(args) > 3:
            for arg in args:
                compounds = pcp.get_compounds(arg, 'name')
                if compounds:
                    bytes = UserCog.monograph(arg)
                    file = discord.File(fp=bytes, filename=f'Molecule.png')
                else:
                    mol = Chem.MolFromSmiles(arg)
                    if mol:
                        try:
                            bytes = UserCog.to_smiles(arg)
                        except Exception as e:
                            await ctx.send(e)
                        file = discord.File(fp=bytes, filename=f'Molecule.png')
                    else:
                        results = search(arg, 'web')[0]['link'], inline=False
                        embed = discord.Embed()
                        for result in results:
                            embed.add_field(name=result['title'], value=result['link'], inline=False)
                        await ctx.send('Invalid string. Enter a molecule name or SMILES string.')
                await ctx.send(file=file)
        if len(args) == 2:
            bytes, similarity = UserCog.monograph_two(args)
            file = discord.File(fp=bytes, filename=f'Molecule.png')
            await ctx.send(file=file)
        if len(args) == 3:
            bytes, similarity = UserCog.monograph_two([args[0], args[1]])
            refbytes, refsimilarity = UserCog.monograph_three([args[2], args[1]])
            img_bytes = Image.open(bytes)
            ref_img_bytes = Image.open(refbytes)
            width, height = img_bytes.size
            new_image = Image.new('RGBA', (width, 2 * height), (0, 0, 0, 0))
            new_image.paste(img_bytes, (0, 0), img_bytes)
            new_image.paste(ref_img_bytes, (0, height), ref_img_bytes)
            draw = ImageDraw.Draw(new_image)
            difference = abs(similarity - refsimilarity)
            modified_image = UserCog.photoshop(draw, new_image, f'-{(difference * 100):.2f}%')
            file = discord.File(fp=modified_image, filename=f'Molecule.png')
            await ctx.send(file=file)


    @commands.command(name='info')
    async def info(self, ctx: commands.Context, argument: str = None):
        if argument is None:
            await ctx.send('Please provide the type of information you need. Options are: `emoji`, `msds`.')
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

    @commands.command(name='smiles')
    async def smiles(self, ctx: commands.Context, *, chemical_name: str):
        compounds = pcp.get_compounds(chemical_name, 'name')
        if not compounds:
            await ctx.send(f'No compound found for the name {chemical_name}')
            return
        compound = compounds[0]
        canonical_smiles = compound.canonical_smiles
        await ctx.send(f'The canonical SMILES for {chemical_name} is: {canonical_smiles}')

    @commands.command(name='get')
    async def get(self, ctx, *, query: str):
        try:
            results = search(query, 'web')
        except Exception as e:
            await ctx.send(e)
        embed = discord.Embed()
        for result in results:
            embed.add_field(name=result['title'], value=result['link'], inline=False)
        await ctx.send(embed=embed)

    @commands.Cog.listener()
    async def on_message(self, message):
        if message.author.bot:
            return
        prefs = self.user_translation_preferences.get(message.author.id)
        if prefs:
            target_lang_code, source_lang_code = prefs
            try:
                translated = self.translator.translate(message.content, dest=target_lang_code, src=source_lang_code)
                await message.channel.send(f'Translated ({source_lang_code} -> {target_lang_code}): {translated.text}')
            except Exception as e:
                await message.channel.send(f'Error translating message: {e}')

    def get_language_code(self, language_name):
        # Convert language name to language code
        language_name = language_name.lower()
        for lang_code, lang_name in LANGUAGES.items():
            if lang_name.lower() == language_name:
                return lang_code
        return None

    @commands.command()
    async def languages(self, ctx):
        # Create a list of supported languages
        supported_languages = ', '.join(LANGUAGES.values())
        await ctx.send(f'Supported languages are:\n{supported_languages}')

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

    async def fetch_image_src(self, query):
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
                    image1 = Image.open(io.BytesIO(image_bytes))
                    width, height = image1.size
                    draw = ImageDraw.Draw(image1)
                    font_size = min(width, height) // 10  # starting font size (adjust as needed)
                    font = ImageFont.load_default(font_size)
                    text_bbox = draw.textbbox((0, 0), query, font=font)
                    text_width = text_bbox[2] - text_bbox[0]
                    text_height = text_bbox[3] - text_bbox[1]
                    while text_width > width - 20:
                        font_size -= 1
                        font = ImageFont.load_default(font_size)
                        text_bbox = draw.textbbox((0, 0), query, font=font)
                        text_width = text_bbox[2] - text_bbox[0]
                        text_height = text_bbox[3] - text_bbox[1]
                    text_x = ((width - text_width) / 2)
                    text_y = height - text_height - 60
                    draw.text((text_x, text_y), query, font=font, fill=(152,251,203))
                    draw.text((text_x + 2, text_y - 2), query, font=font, fill=(152,251,203))
                    draw.text((text_x + 4, text_y - 4), query, font=font, fill=(152,251,203))
                    output_buffer = io.BytesIO()
                    image1.save(output_buffer, format='PNG')
                    output_buffer.seek(0)
                    file = discord.File(output_buffer, filename="image.png")
                    return file
                else:
                    return "No src attribute found in the <img> element"
            else:
                return "No <img> element found with the specified CSS path"
        finally:
            driver.quit()

    @commands.command(name='search')
    async def search(self, ctx, *, query: str):
        # Ensure the event loop is running
        loop = asyncio.get_event_loop()
        try:
            #img_src = await loop.run_in_executor(None, lambda: asyncio.run(UserCog.fetch_image_src(self, query)))
            try:
                img_src = await loop.run_in_executor(None, lambda: asyncio.run(UserCog.fetch_wm_image_src(self, query)))
            except Exception as e:
                await ctx.send(e)

            await ctx.send(file=img_src)
        except Exception as e:
            await ctx.send(f'{query} is an unknown molecule.')

    def fetch_wm_image_src(self, query):
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
                    image1 = Image.open(io.BytesIO(image_bytes))
                    width, height = image1.size
                    draw = ImageDraw.Draw(image1)
                    diagonal = math.sqrt(width**2 + height**2)
                    font_size = int(diagonal / 15)  # Adjust this factor as needed
                    try:
                        font = ImageFont.truetype("arial.ttf", font_size)
                    except IOError:
                        font = ImageFont.load_default()
                    bbox = draw.textbbox((0, 0), query, font=font)
                    text_width = bbox[2] - bbox[0]
                    text_height = bbox[3] - bbox[1]
                    text_x = (width - text_width) / 2
                    text_y = (height - text_height) / 2
                    watermark_image = Image.new("RGBA", image.size, (0, 0, 0, 0))
                    watermark_draw = ImageDraw.Draw(watermark_image)
                    watermark_draw.text((text_x, text_y), query, font=font, fill=(255, 255, 255, 64))
                    watermark_image = watermark_image.rotate(45, expand=True)
                    wm_width, wm_height = watermark_image.size
                    paste_x = (width - wm_width) / 2
                    paste_y = (height - wm_height) / 2
                    image1.paste(watermark_image, (int(paste_x), int(paste_y)), watermark_image)
                    output_buffer = io.BytesIO()
                    image1.save(output_buffer, format='PNG')
                    output_buffer.seek(0)
                    file = discord.File(output_buffer, filename="image.png")
                    return file
                else:
                    return "No src attribute found in the <img> element"
            else:
                return "No <img> element found with the specified CSS path"
        finally:
            driver.quit()

async def setup(bot: commands.Bot):
    await bot.add_cog(UserCog(bot))
