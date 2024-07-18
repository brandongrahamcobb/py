
# cogs/chemistry.py

from discord.ext import commands
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import SimilarityMaps
from rdkit.DataStructs import FingerprintSimilarity, TanimotoSimilarity

import discord
import io
import pubchempy as pcp

class Chemistry(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self._last_member = None

    @commands.command(description='Compare molecules, fetch molecules.')
    async def draw(self, ctx: commands.Context, *args) -> None:
        if len(args) == 1:
             mol = Chem.MolFromSmiles(args[0])
             img = Draw.MolToImage(mol, size=(512, 512))
             mol = io.BytesIO()
             img.save(mol, format='PNG')
             mol.seek(0)
             file = discord.File(fp=mol, filename=f'Molecule.png')
             await ctx.send(file=file)
        compounds = pcp.get_compounds(args[0], 'name')
        compound_data = compounds[0].to_dict(properties=['isomeric_smiles'])
        label_1 = compounds[0].synonyms[0]
        mol = Chem.MolFromSmiles(compound_data['isomeric_smiles'])
        other_compounds = pcp.get_compounds(args[1], 'name')
        other_compound_data = other_compounds[0].to_dict(properties=['isomeric_smiles'])
        label_2 = other_compounds[0].synonyms[0]
        refmol = Chem.MolFromSmiles(other_compound_data['isomeric_smiles'])
        d2d = Draw.MolDraw2DCairo(512, 512)
        d2d2 = Draw.MolDraw2DCairo(512, 512)
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(refmol, 2, nBits=2048)
        similarity = TanimotoSimilarity(fp1, fp2)
        fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(refmol, mol, lambda m, i: SimilarityMaps.GetMorganFingerprint(m, i, radius=2, fpType='bv', nBits=4096), draw2d=d2d)
        d2d.FinishDrawing()
        buf = io.BytesIO(d2d.GetDrawingText())
        image1 = Image.open(buf)
        image1 = image1.convert("RGBA")
        fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(mol, refmol, lambda m, i: SimilarityMaps.GetMorganFingerprint(m, i, radius=2, fpType='bv', nBits=4096), draw2d=d2d2)
        d2d2.FinishDrawing()
        buf = io.BytesIO(d2d2.GetDrawingText())
        image2 = Image.open(buf)
        image2 = image2.convert("RGBA")
        width, height = image2.size
        new_image = Image.new('RGBA', (width + width, height), (0, 0, 0, 0))
        new_image.paste(image1, (0, 0), image1)
        new_image.paste(image2, (width, 0), image2)
        draw = ImageDraw.Draw(new_image)
        for x in range(new_image.width):
            for y in range(new_image.height):
                pixel = new_image.getpixel((x, y))
                if sum(pixel) >= 1020:
                    draw.point((x, y), fill=(0, 0, 0, 0))
        font = ImageFont.load_default(40)
        draw.text((((width - draw.textlength(label_1, font=font)) / 2) + 4, 14), label_1, fill="gray", font=font)
        draw.text(((width + (width - draw.textlength(label_2, font=font)) / 2) + 4, 14), label_2, fill="gray", font=font)
        draw.text(((width - draw.textlength(label_1, font=font)) / 2, 10), label_1, fill="black", font=font)
        draw.text((width + (width - draw.textlength(label_2, font=font)) / 2, 10), label_2, fill="black", font=font)
        text = f"Tanimoto Similarity: {similarity:.5f}"
        draw.text((((new_image.width - draw.textlength(text, font=font)) / 2) + 4, height - 56), text, fill="gray", font=font)
        draw.text(((new_image.width - draw.textlength(text, font=font)) / 2, height - 60), text, fill="black", font=font)
        output_buffer = io.BytesIO()
        new_image.save(output_buffer, format='PNG')
        output_buffer.seek(0)
        file = discord.File(fp=output_buffer, filename=f'Molecule.png')
        await ctx.send(file=file)

async def setup(bot):
    await bot.add_cog(Chemistry(bot))
