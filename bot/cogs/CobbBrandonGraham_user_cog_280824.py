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
from discord import app_commands
from discord.ext import commands
from typing import List

import bot.utils.helpers as lucy

import discord
import io
import os
import pubchempy as pcp
import rdkit
import requests
import shlex
import traceback

class UserCog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.bot.config = bot.config
        self.stacks = {}

    @commands.Cog.listener()
    async def on_message(self, message: discord.Message):
        if message.author.bot: # or message.channel.id == '962009752488013834':
            return

    def get_user_stack(self, user_id: int):
        return self.stacks.get(user_id, [])

    def set_user_stack(self, user_id: int, molecule_names: List[str]):
        mol_objects = [lucy.get_mol(name) for name in molecule_names]
        self.stacks[user_id] = mol_objects

    @app_commands.command(name='draw', description='Draw a molecule or compare molecules by their names.')
#    @app_commands.describe(molecules='')
    async def draw(self, interaction: discord.Interaction, molecules: str) -> None:
        try:
            await interaction.response.defer(ephemeral=True)
            args = shlex.split(molecules)
            if len(args) == 1:
                mol = lucy.get_mol(args[0])
                if mol is None:
                    embed = 'Invalid molecule name or structure.'
                    await interaction.followup.send(embed=embed)
                    return
                image = lucy.draw_watermarked_molecule(mol)
                
                await interaction.followup.send(file=discord.File(image, f'{args[0]}.png'))
                return
            pairs = lucy.unique_pairs(args)
            if not pairs:
                embed = 'No valid pairs found.'
                await interaction.followup.send(embed=embed)
                return
            for pair in pairs:
                mol = lucy.get_mol(pair[0])
                refmol = lucy.get_mol(pair[1])
                if mol is None or refmol is None:
                    embed = f'One or both of the molecules {pair[0]} or {pair[1]} are invalid.'
                    await interaction.followup.send(embed=embed)
                    continue
                #image = lucy.draw_watermarked_molecule(mol)
                #proximity = lucy.get_proximity(default=mol, input=refmol)
                #embed = f'The Tanimoto Similarity of {lucy.get_molecule_name(mol)} and {lucy.get_molecule_name(refmol)} is {proximity:.2f}'
                #await interaction.followup.send(embed=embed)
                fingerprints = [
                   lucy.draw_fingerprint([mol, refmol]),
                    lucy.draw_fingerprint([refmol, mol])
                ]
                combined_image = lucy.combine(fingerprints[0], fingerprints[1])
                await interaction.followup.send(file=discord.File(combined_image, f'molecule_comparison.png'))
        except Exception as e:
            await interaction.followup.send(f'An error occurred: {str(e)}')

    @app_commands.command(name='get', hidden=True)
    @commands.is_owner()
    async def get(self, interaction: discord.Interaction, image: str):
        result = get_images(image)
        await interaction.response.send_message(result)

    @app_commands.command(name='info', hidden=True)
    async def info(self, interaction: discord.Interaction, argument):
       if argument is None:
           await interaction.response.send_message('The command is used !info <MOLECULE> or !info :emoji:')
           return
       if isinstance(argument, discord.Emoji):
            try:
               await interaction.response.send_message(lucy.get_emoji(argument))
            except Exception as e:
                await interaction.response.send_message(e)
       else:
            try:
               embed = lucy.get_sds(argument)
               await interaction.response.send_message(embed=embed)
            except Exception as e:
               await interaction.response.send_message(e)

    @app_commands.command(name='logp', hidden=True)
    async def logp(self, interacton: discord.Interaction, molecules: str):
        try:
            args = shlex.split(molecules)
            for arg in args:
                compounds = pcp.get_compounds(arg, 'name')
                compound_data = compounds[0].to_dict(properties=['isomeric_smiles'])
                mol = Chem.MolFromSmiles(compound_data['isomeric_smiles'])
                log_p = Crippen.MolLogP(mol)
                await interaction.response.send_message(f'Your octanol:water coefficient is: {log_p}')
        except Exception as e:
            await interaction.response.send_message(f'Error fetching data: {e}')

    @app_commands.command(name='search', hidden=True)
    async def search(self, interaction: discord.Interaction *, arg=None):
        try:
            async with ctx.typing():
                try:
                    watermarked_image = lucy.gsrs(arg)
                    with io.BytesIO() as image_binary:
                        watermarked_image.save(image_binary, format='PNG')
                        image_binary.seek(0)
                        await interaction.response.send_message(file=discord.File(fp=image_binary, filename='watermarked_image.png'))
                except:
                    await interaction.response.send_message(traceback.format_exc())
        except:
            if arg is None:
                return
            await interaction.response.send_message(f'{arg} is an unknown molecule.')

    @app_commands.command(name='smiles', hidden=True)
    async def smiles(self, ctx: commands.Context, *, chemical_name: str):
        async with ctx.typing():
            compounds = pcp.get_compounds(chemical_name, 'name')
            compound = compounds[0]
            isomeric_smiles = compound.isomeric_smiles
            await interaction.response.send_message(f'The isomeric SMILES for {chemical_name} is: {isomeric_smiles}')

    @app_commands.command(name='stack', description='Submit one molecule to test or multiple to set a stack.')
    async def stack(self, ctx: commands.Context, *molecules):
        async with ctx.typing():
            if not molecules:
                await interaction.response.send_message([lucy.get_molecule_name(mol) for mol in self.get_user_stack(ctx.author.id)])
            if ctx.author.id in self.stacks and len(molecules) == 1:
                defender = max(self.get_user_stack(ctx.author.id), key=lambda mol: lucy.get_proximity(mol, lucy.get_mol(molecules[0])))
                proximity = lucy.get_proximity(lucy.get_mol(molecules[0]), defender)
                await interaction.response.send_message(
                    f'{molecules[0]} | {proximity:.3f}% | {lucy.get_molecule_name(defender)}\n'
                )
            elif len(molecules) > 1:
                self.set_user_stack(user_id=ctx.author.id, molecule_names=molecules)
                await interaction.response.send_message('Stack overwritten.')

async def setup(bot: commands.Bot):
    await bot.add_cog(UserCog(bot))
