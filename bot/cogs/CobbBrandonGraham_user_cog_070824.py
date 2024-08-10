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

from bot.utils.helpers import google_search
from bot.utils.helpers import get_emoji_info
from bot.utils.helpers import get_msds_info
from bot.utils.helpers import get_version
from bot.utils.helpers import load_config
from bot.utils.helpers import same_hashes
from discord.ext import commands
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import SimilarityMaps
from rdkit.DataStructs import FingerprintSimilarity, TanimotoSimilarity

import discord
import json
import io
import pubchempy as pcp

class UserCog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        config = load_config()
        self.google_json_key = config.get('google_json_key')
        self.my_cse_id = '4231a573aa34241e4'
        self.your_cse_id = '25fb95395119b40f0'
        self.new_cse_id = 'e1753e22095c74f10'

    def embedded_list(item: str):
        embed = discord.Embed()
        results = MyCog.google_search(item, self.google_json_key, self.my_cse_id, num=10)
        for result in results:
            embed.add_field(name=result['title'], value=result['link'], inline=False)
        return embed

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
        image1 = image1.convert("RGBA")
        width, height = image1.size
        draw = ImageDraw.Draw(image1)
        for x in range(image1.width):
            for y in range(image1.height):
                pixel = image1.getpixel((x, y))
                if sum(pixel) >= 1020:
                    draw.point((x, y), fill=(0, 0, 0, 0))
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
        draw.text((text_x + 4, text_y + 4), text, font=font, fill=(152,251,203))
        output_buffer = io.BytesIO()
        image1.save(output_buffer, format='PNG')
        output_buffer.seek(0)
        return output_buffer

    def smiles(molecule: str):
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
        return img_bytes

    @commands.command(description='Compare molecules, fetch molecules.')
    async def draw(self, ctx: commands.Context, *args) -> None:
        if len(args) == 1:
            compounds = pcp.get_compounds(args[0], 'name')
            if compounds:
                bytes = UserCog.monograph(args[0])
                file = discord.File(fp=bytes, filename=f'Molecule.png')
            else:
                mol = Draw.MolFromSmiles(args[0])
                if mol:
                    bytes = UserCog.smiles(args[0])
                    file = discord.File(fp=bytes, filename=f'Molecule.png')
                    await ctx.send(file=file)
                else:
                    embed = self.embedded_list(args[0])
                    message = await ctx.send(embed=embed)
                    if not message:
                        await ctx.send('Invalid string. Enter a molecule name or SMILES string.')
            await ctx.send(file=file)
        if len(args) == 2:
            d2d1 = Draw.MolDraw2DCairo(512, 512)
            d2d2 = Draw.MolDraw2DCairo(512, 512)
            compounds = pcp.get_compounds(args[0], 'name')
            if compounds:
                label_1 = compounds[0].synonyms[0]
                compound_data = compounds[0].to_dict(properties=['isomeric_smiles'])
                mol = Chem.MolFromSmiles(compound_data['isomeric_smiles'])
            else:
                label_1 = 'Unknown'
                mol = Draw.MolFromSmiles(args[0])
            ref_compounds = pcp.get_compounds(args[1], 'name')
            if ref_compounds:
                label_2 = ref_compounds[0].synonyms[0]
                ref_compound_data = ref_compounds[0].to_dict(properties=['isomeric_smiles'])
                refmol = Chem.MolFromSmiles(ref_compound_data['isomeric_smiles'])
            else:
                label_2 = 'Unknown'
                refmol = Draw.MolFromSmiles(args[1])
                await ctx.send('Invalid first string. Enter a molecule name or SMILES string.')
                await ctx.send('Invalid second string. Enter a molecule name or SMILES string.')
            fp1 = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            fp2 = AllChem.GetMorganFingerprintAsBitVect(refmol, 2, nBits=2048)
            similarity = TanimotoSimilarity(fp1, fp2)
            fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(refmol, mol, lambda m, i: SimilarityMaps.GetMorganFingerprint(m, i, radius=2, fpType='bv', nBits=8192), draw2d=d2d1)
            d2d1.FinishDrawing()
            buf = io.BytesIO(d2d1.GetDrawingText())
            image1 = Image.open(buf)
            image1 = image1.convert("RGBA")
            fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(mol, refmol, lambda m, i: SimilarityMaps.GetMorganFingerprint(m, i, radius=2, fpType='bv', nBits=8192), draw2d=d2d2)
            d2d2.FinishDrawing()
            buf = io.BytesIO(d2d2.GetDrawingText())
            image2 = Image.open(buf)
            image2 = image2.convert("RGBA")
            width, height = image2.size
            new_image = Image.new('RGBA', (width + width, height), (0, 0, 0, 0))
            new_image.paste(image1, (0, 0), image1)
            new_image.paste(image2, (width, 0), image2)
            original = new_image
            transparent_version = original.copy()
            alpha = transparent_version.split()[3]
            alpha = alpha.point(lambda p: min(p, 128))
            transparent_version.putalpha(alpha)
            combined = Image.new("RGBA", original.size)
            combined.paste(transparent_version, (10, 10), mask=transparent_version)
            combined.paste(original, (0, 0))
            draw = ImageDraw.Draw(combined)
            font = ImageFont.load_default(40)
            for x in range(combined.width):
                for y in range(combined.height):
                    pixel = combined.getpixel((x, y))
                    if sum(pixel) >= 1020:
                        draw.point((x, y), fill=(0, 0, 0, 0))
            draw.text(((width - draw.textlength(label_1, font=font)) / 2, 10), label_1, fill=(152,251,203), font=font)
            draw.text((width + (width - draw.textlength(label_2, font=font)) / 2, 10), label_2, fill=(152,251,203), font=font)
            draw.text((((width - draw.textlength(label_1, font=font)) / 2) + 4, 6), label_1, fill=(152,251,203), font=font)
            draw.text(((width + (width - draw.textlength(label_2, font=font)) / 2) + 4, 6), label_2, fill=(152,251,203), font=font)
            draw.text(((0.5 * width) + (width - draw.textlength(f'LUCY v{get_version()}', font=font)) / 2, 10), f'LUCY v{get_version()}', fill="yellow", font=font)
            text = f"Tanimoto Similarity: {similarity:.5f}"
            draw.text(((combined.width - draw.textlength(text, font=font)) / 2, height - 60), text, fill="green", font=font)
            output_buffer = io.BytesIO()
            combined.save(output_buffer, format='PNG')
            output_buffer.seek(0)
            file = discord.File(fp=output_buffer, filename=f'Molecule.png')
            await ctx.send(file=file)

    @commands.command(name='info')
    async def info(self, ctx: commands.Context, info_type: str = None, *, argument: str = None):
        if info_type is None:
            await ctx.send("Please provide the type of information you need. Options are: `emoji`, `msds`.")
            return
        if info_type.lower() == "emoji":
            if argument is None:
                await ctx.send("Please provide the Unicode emoji character.")
            else:
                try:
                    await get_emoji_info(ctx, argument)
                except Exception as e:
                    await ctx.send(e)
        elif info_type.lower() == "msds":
            if argument is None:
                await ctx.send("Please provide the molecule name.")
            else:
                try:
                    await get_msds_info(ctx, argument)
                except Exception as e:
                    await ctx.send(e)
        else:
            await ctx.send("Invalid info type. Options are: `emoji`, `msds`.")

    @commands.command(name='smiles')
    async def smiles(self, ctx: commands.Context, *, chemical_name: str):
        compounds = pcp.get_compounds(chemical_name, 'name')
        if not compounds:
            await ctx.send(f"No compound found for the name '{chemical_name}'")
            return
        compound = compounds[0]
        canonical_smiles = compound.canonical_smiles
        await ctx.send(f"The canonical SMILES for '{chemical_name}' is: {canonical_smiles}")

#    @commands.command(name='get')
 #   async def get(self, ctx, *, query: str):
  #      results = MyCog.google_search(f'{query}', self.google_json_key, self.new_cse_id, num=4)
   #     embed = discord.Embed()
    #    for result in results:
     #       embed.add_field(name=result['title'], value=result['link'], inline=False)
      #  await ctx.send(embed=embed)

async def setup(bot: commands.Bot):
    await bot.add_cog(UserCog(bot))
