''' my_cog.py
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
from discord.ext import commands
from PIL import Image, ImageFont, ImageDraw
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Draw
from rdkit.Chem.Draw import rdMolDraw2D, SimilarityMaps
from rdkit.DataStructs import FingerprintSimilarity, TanimotoSimilarity
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)
from rdkit.Chem.Draw import IPythonConsole
import discord
from io import BytesIO
import io
import math
import PIL
import pubchempy as pcp
import rdkit
import traceback

class MyCog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot

    @commands.Cog.listener()
    async def on_command_error(self, ctx: commands.Context, error):
        if isinstance(error, commands.CommandOnCooldown):
            await ctx.send(f'Command is on cooldown. Try again in {error.retry_after:.2f} seconds.')
        else:
            raise error

    def add_watermark(self, image: BytesIO, watermark_text: str) -> BytesIO:
        RGB_image = Image.open(image)
        RGBA_image = RGB_image.convert('RGBA')
        draw = ImageDraw.Draw(RGBA_image)
        width, height = RGBA_image.size
        diagonal = math.sqrt(width**2 + height**2)
        font_size = int(diagonal / 15)  # Adjust this factor as needed
        font = ImageFont.load_default(font_size)
        bbox = draw.textbbox((0, 0), watermark_text, font=font)
        text_width = bbox[2] - bbox[0]
        text_height = bbox[3] - bbox[1]
        text_x = (width - text_width) / 2
        text_y = (height - text_height) / 2
        watermark_image = Image.new('RGBA', RGBA_image.size, (0, 0, 0, 0))
        watermark_draw = ImageDraw.Draw(watermark_image)
        watermark_draw.text((text_x, text_y), watermark_text, font=font, fill=(255, 255, 255, 64))
        mask = watermark_image.split()[3]
        RGBA_image.paste(watermark_image, (0, 0), mask)
        output = BytesIO()
        RGBA_image.save(output, format= 'PNG')
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
        output = BytesIO()
        combined_img.save(output, format='PNG')
        output.seek(0)
        return output

    def draw_watermarked_fingerprint(self, queries) -> BytesIO:
        molecules = []
        for smiles in queries:
            mol = Chem.MolFromSmiles(self.get_smiles(smiles))
            molecules.append(mol)
        compounds = pcp.get_compounds(self.get_smiles(queries[0]), 'smiles')
        resolved_name = compounds[0].synonyms[0]
        d2d = rdMolDraw2D.MolDraw2DCairo(512, 512)
        d2d.prepareMolsBeforeDrawing = False
        Options = d2d.drawOptions()
        Options.prepareMolsBeforeDrawing = False
        Options.includeMetadata = False
        d2d.SetDrawOptions(Options)
        mol1 = rdMolDraw2D.PrepareMolForDrawing(molecules[0], kekulize=True)
        mol1.UpdatePropertyCache(False)
        mol2 = rdMolDraw2D.PrepareMolForDrawing(molecules[1], kekulize=True)
        mol2.UpdatePropertyCache(False)
        fp1 = AllChem.GetMorganFingerprintAsBitVect(molecules[0], 2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(molecules[1], 2, nBits=2048)
        similarity = TanimotoSimilarity(fp1, fp2)
        fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(mol1, mol2, lambda m, i: SimilarityMaps.GetMorganFingerprint(m, i, radius=2, fpType='bv', nBits=8192), draw2d=d2d, drawingOptions=Options)
        d2d.FinishDrawing()
        drawing = d2d.GetDrawingText()
        output = self.add_watermark(BytesIO(drawing), watermark_text=resolved_name)
        return output

    def draw_watermarked_molecule(self, query: str) -> BytesIO:
        smiles = self.get_smiles(query)
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                compounds = pcp.get_compounds(smiles, 'smiles')
                resolved_name = compounds[0].synonyms[0]
                d2d = rdMolDraw2D.MolDraw2DCairo(512, 512)
                Options = d2d.drawOptions()
                d2d.SetDrawOptions(Options)
                rdMolDraw2D.SetDarkMode(Options)
                mol1 = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=True)
                mol1.UpdatePropertyCache(False)
                d2d.DrawMolecule(mol1)
                d2d.FinishDrawing()
                drawing = d2d.GetDrawingText()
                output = self.add_watermark(BytesIO(drawing), watermark_text=resolved_name)
                return output
            else:
                return None
        else:
            return None

    def get_smiles(self, query: str) -> str:
        mol = Chem.MolFromSmiles(query) #or Chem.MolFromInchl(query)
        if mol:
            return query
        else:
            results = pcp.get_compounds(query, 'name')
            if results:
                return results[0].isomeric_smiles
            else:
                return None

    @commands.command(name='draw')
    async def draw(self, ctx, *queries) -> None:
        cooldown_manager = self.bot.get_cog("AdminCog")
        if cooldown_manager:
            cooldown = cooldown_manager.get_cooldown(ctx.channel.id, ctx.command.name)
            if cooldown > 0:
                bucket = ctx.command._buckets.get_bucket(ctx)
                retry_after = bucket.update_rate_limit()
                if retry_after:
                    await ctx.send(f"Command is on cooldown. Try again in {retry_after:.2f} seconds.")
                    return
        try:
            molecules = []
            for query in queries:
                image = self.draw_watermarked_molecule(query)
                if image:
                    molecules.append(image)
            if len(molecules) == 1 or len(molecules) > 3:
                for i, image in enumerate(molecules):
                    await ctx.send(file=discord.File(image, f'molecule_{i+1}.png'))
            elif len(molecules) == 2:
                fingerprints = []
                fingerprints.append(self.draw_watermarked_fingerprint(list(reversed(queries))))
                fingerprints.append(self.draw_watermarked_fingerprint(queries))
                combined_image = self.combine(fingerprints[0], fingerprints[1])
                await ctx.send(file=discord.File(combined_image, f'molecule_comparison.png'))
            else:
                await ctx.send('No valid molecules found.')
        except Exception as e:
            await ctx.send(f'{e} {traceback.print_exc()}')

    @commands.command(name='logp')
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

async def setup(bot: commands.Bot):
    await bot.add_cog(MyCog(bot))
