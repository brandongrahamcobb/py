''' info_cog.py The purpose of this program is to provide functionality to a caffeine replacement therapy bot on Discord.
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

class Molecular(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config

    @commands.hybrid_command(name='draw', description='Perform various operations on molecules')
    async def molecule(self, ctx: commands.Context, option=None, *, molecules=None, quantity: int = 1):
        try:
            if ctx.interaction:
                await ctx.interaction.response.defer(ephemeral=True)
            if option == 'compare':
                if not molecules:
                    await ctx.send('No molecules provided.')
                    return
                args = shlex.split(molecules)
                pairs = helpers.unique_pairs(args)
                if not pairs:
                    embed = discord.Embed(description='No valid pairs found.')
                    await ctx.send(embed=embed)
                    return
                for pair in pairs:
                    mol = helpers.get_mol(pair[0])
                    refmol = helpers.get_mol(pair[1])
                    if mol is None or refmol is None:
                        embed = discord.Embed(description=f'One or both of the molecules {pair[0]} or {pair[1]} are invalid.')
                        await ctx.send(embed=embed)
                        continue
                    fingerprints = [
                        helpers.draw_fingerprint([mol, refmol]),
                        helpers.draw_fingerprint([refmol, mol])
                    ]
                    combined_image = helpers.combine(fingerprints, reversed(pair))
                    await ctx.send(file=discord.File(combined_image, f'molecule_comparison.png'))
            elif option == 'glow':
                if not molecules:
                    await ctx.send('No molecules provided.')
                    return
                args = shlex.split(molecules)
                fingerprints = []
                names = []
                molecule = helpers.get_mol(args[0])
                for _ in range(quantity):
                    names.append(args[0])
                    fingerprints.append(helpers.draw_fingerprint([molecule, molecule]))
                combined_image = helpers.combine(fingerprints, names)
                await ctx.send(file=discord.File(combined_image, f'molecule_comparison.png'))
            elif option == 'gsrs':
                if not molecules:
                    await ctx.send('No molecules provided.')
                    return
                args = shlex.split(molecules)
                for molecule_name in args:
                    if molecule_name is None:
                        await ctx.send(f'{molecule_name} is an unknown molecule.')
                        continue
                    watermarked_image = helpers.gsrs(molecule_name)
                    with io.BytesIO() as image_binary:
                        watermarked_image.save(image_binary, format='PNG')
                        image_binary.seek(0)
                        await ctx.send(file=discord.File(fp=image_binary, filename='watermarked_image.png'))
            elif option == 'shadow':
                if not molecules:
                    await ctx.send('No molecules provided.')
                    return
                args = shlex.split(molecules)
                mol = helpers.get_mol(args[0])
                if mol is None:
                    embed = discord.Embed(description='Invalid molecule name or structure.')
                    await ctx.send(embed=embed)
                    return
                image = helpers.draw_watermarked_molecule(mol)
                await ctx.send(file=discord.File(image, f'{args[0]}.png'))
            else:
                await ctx.send('Invalid option. Use `draw`, `compare`, `draw_alt`, or `draw_alt_2`.')
        except Exception as e:
            print(f'An error occurred: {traceback.format_exc()}')

async def setup(bot: commands.Bot):
    await bot.add_cog(Molecular(bot))
