''' chem_cog.py The purpose of this program is to provide core chemistry functionality to a bot on Discord.
    Copyright (C) 2024  github.com/brandongrahamcobb

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

from discord import app_commands, Embed
from discord.ext import commands

import bot.utils.helpers as helpers

import discord
import io
import shlex
import traceback

class ChemCog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config

    @commands.hybrid_command(name='draw', description='Draw a molecule or compare molecules by their names.')
    async def draw(self, ctx: commands.Context, *, molecules: str) -> None:
        try:
            if ctx.interaction:
                await ctx.interaction.response.defer(ephemeral=True)
            args = shlex.split(molecules)
            if len(args) == 1:
                mol = helpers.get_mol(args[0])
                if mol is None:
                    embed = 'Invalid molecule name or structure.'
                    await ctx.send(embed=embed)
                    return
                image = helpers.draw_watermarked_molecule(mol)
                await ctx.send(file=discord.File(image, f'{args[0]}.png'))
                return
            pairs = helpers.unique_pairs(args)
            if not pairs:
                embed = 'No valid pairs found.'
                await ctx.send(embed=embed)
                return
            for pair in pairs:
                mol = helpers.get_mol(pair[0])
                refmol = helpers.get_mol(pair[1])
                if mol is None or refmol is None:
                    embed = f'One or both of the molecules {pair[0]} or {pair[1]} are invalid.'
                    await ctx.send(embed=embed)
                    continue
                #image = helpers.draw_watermarked_molecule(mol)
                #proximity = helpers.get_proximity(default=mol, input=refmol)
                #embed = f'The Tanimoto Similarity of {helpers.get_molecule_name(mol)} and {helpers.get_molecule_name(refmol)} is {proximity:.2f}'
                #await interaction.followup.send(embed=embed)
                fingerprints = [
                    helpers.draw_fingerprint([mol, refmol]),
                    helpers.draw_fingerprint([refmol, mol])
                ]
                combined_image = helpers.combine(fingerprints, reversed(pair))
                await ctx.send(file=discord.File(combined_image, f'molecule_comparison.png'))
        except Exception as e:
            await ctx.send(f'An error occurred: {str(e)}')

    @commands.hybrid_command(name='altdraw', description='Draws chemicals with transparency.')
    async def draw_alt(self, ctx: commands.Context, *, molecules: str):
        try:
            if ctx.interaction:
                 await ctx.interaction.response.defer(ephemeral=True)
            args = shlex.split(molecules)
            try:
                for arg in args:
                    watermarked_image = helpers.gsrs(molecules)
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

    @commands.hybrid_command(name='altdraw2', description='Draws chemicals with glow.')
    async def draw_alt_alt(self, ctx: commands.Context, molecules: str, quantity: int):
        fingerprints = []
        names = []
        if not molecules or molecules == 'cocaine' or 'cocaethylene' or 'coke' and quantity > 100:
            if ctx.interaction:
                await ctx.interaction.response.defer(ephemeral=True)
            args = shlex.split(molecules)
            molecule = helpers.get_mol(args[0])
            for mol in range(quantity):
                names.append(args[0])
                fingerprints.append(helpers.draw_fingerprint([molecule, molecule]))
            combined_image = helpers.combine(fingerprints, names)
            await ctx.send(file=discord.File(combined_image, f'molecule_comparison.png'))

async def setup(bot: commands.Bot):
    await bot.add_cog(ChemCog(bot))
