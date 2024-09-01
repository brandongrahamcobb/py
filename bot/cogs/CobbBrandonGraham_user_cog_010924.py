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
from rdkit import Chem
from rdkit.Chem import Crippen
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
        self.lysergic_acid_diethylamide = lucy.get_mol('lysergic acid diethylamide')
        self.quetiapine = lucy.get_mol('quetiapine')

    @commands.Cog.listener()
    async def on_message(self, message: discord.Message):
        if message.author.bot: # or message.channel.id == '962009752488013834':
            return

    def get_user_stack(self, user_id: int):
        return self.stacks.get(user_id, [])

    def set_user_stack(self, user_id: int, molecule_names: List[str]):
        mol_objects = [lucy.get_mol(name) for name in molecule_names]
        self.stacks[user_id] = mol_objects

    @commands.command(name='analog', description='Submit one molecule to test or multiple to set a stack.')
    async def analog(self, ctx: commands.Context, *args):
        if args:
            try:
                self.set_user_stack(user_id=ctx.author.id, molecule_names=args)
                await ctx.send([lucy.get_molecule_name(mol) for mol in self.get_user_stack(ctx.author.id)])
                await ctx.send('Stack overwritten. Which molecule would you like to know about?')
                while True:
                    response = await self.bot.wait_for(
                        'message',
                        timeout=600.0,
                         check=lambda message: message.author == ctx.author and message.channel == ctx.channel
                    )
                    try:
                        molecule = lucy.get_mol(response.content)
                        defender = max(self.get_user_stack(ctx.author.id), key=lambda mol: lucy.get_proximity(mol, molecule))
                        proximity = lucy.get_proximity(molecule, defender)
                        await ctx.send(f'{lucy.get_molecule_name(molecule)} is most analogous to {lucy.get_molecule_name(defender)} by a degree of {(100 * proximity):.3f}%')
                        if ctx.author.id == 154749533429956608:
                            max_lsd = max(self.get_user_stack(ctx.author.id), key=lambda mol: lucy.get_proximity(mol, self.lysergic_acid_diethylamide))
                            max_q = max(self.get_user_stack(ctx.author.id), key=lambda mol: lucy.get_proximity(mol, self.quetiapine))
                            max_lsd_proximity = lucy.get_proximity(molecule, max_lsd)
                            max_q_proximity = lucy.get_proximity(molecule, max_q)
                            if max_lsd == max_q:
                                await ctx.send(f'{(100 * max_lsd_proximity):.3f}% | L | You\'d have to stop taking {lucy.get_molecule_name(max_lsd)} to take {lucy.get_molecule_name(molecule)} | Q | {(100 * max_q_proximity):.3f}%')
                            await ctx.send(f'{(100 * max_lsd_proximity):.3f}% = {lucy.get_molecule_name(max_lsd)} | L | {lucy.get_molecule_name(molecule)} | Q | {lucy.get_molecule_name(max_q)} = {(100 * max_q_proximity):.3f}%')
                        await ctx.send('Which molecule would you like to know about next?')
                    except:
                        break
            except:
                await ctx.send('Timeout error. Restart the command.')
        else:
             await ctx.send('Invalid option. Submit molecule(s) to build a library.')

    @commands.hybrid_command(name='draw', description='Draw a molecule or compare molecules by their names.')
    async def draw(self, ctx: commands.Context, molecules: str) -> None:
        try:
            if ctx.interaction:
                await ctx.interaction.response.defer(ephemeral=True)
#            args = molecules.split()
            args = shlex.split(molecules)
            if len(args) == 1:
                mol = lucy.get_mol(args[0])
                if mol is None:
                    embed = 'Invalid molecule name or structure.'
                    await ctx.send(embed=embed)
                    return
                image = lucy.draw_watermarked_molecule(mol)
                await ctx.send(file=discord.File(image, f'{args[0]}.png'))
                return
            pairs = lucy.unique_pairs(args)
            if not pairs:
                embed = 'No valid pairs found.'
                await ctx.send(embed=embed)
                return
            for pair in pairs:
                mol = lucy.get_mol(pair[0])
                refmol = lucy.get_mol(pair[1])
                if mol is None or refmol is None:
                    embed = f'One or both of the molecules {pair[0]} or {pair[1]} are invalid.'
                    await ctx.send(embed=embed)
                    continue
                #image = lucy.draw_watermarked_molecule(mol)
                #proximity = lucy.get_proximity(default=mol, input=refmol)
                #embed = f'The Tanimoto Similarity of {lucy.get_molecule_name(mol)} and {lucy.get_molecule_name(refmol)} is {proximity:.2f}'
                #await interaction.followup.send(embed=embed)
                fingerprints = [
                   lucy.draw_fingerprint([mol, refmol]),
                    lucy.draw_fingerprint([refmol, mol])
                ]
                combined_image = lucy.combine(fingerprints[0], pair[1], fingerprints[1], pair[0])
                await ctx.send(file=discord.File(combined_image, f'molecule_comparison.png'))
        except Exception as e:
            await ctx.send(f'An error occurred: {str(e)}')

    @app_commands.command(name='get')
    @commands.is_owner()
    async def get(self, interaction: discord.Interaction, image: str):
        result = get_images(image)
        await interaction.response.send_message(result)

    @app_commands.command(name='emoji')
    async def emoji(self, interaction: discord.Interaction, argument: str):
        try:
           await interaction.response.send_message(lucy.get_emoji(argument))
        except Exception as e:
            await interaction.response.send_message(e)

    @app_commands.command(name='logp')
    async def logp(self, interaction: discord.Interaction, molecules: str):
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

    @app_commands.command(name='msds')
    async def msds(self, interaction: discord.Interaction, argument: str):
        try:
           embed = lucy.get_sds(argument)
           await interaction.response.send_message(embed=embed)
        except Exception as e:
           await interaction.response.send_message(e)

    @app_commands.command(name='search')
    async def search(self, interaction: discord.Interaction, molecules: str):
        try:
            await interaction.response.defer(ephemeral=True)
            args = shlex.split(molecules)
            try:
                for arg in args:
                    watermarked_image = lucy.gsrs(molecules)
                    with io.BytesIO() as image_binary:
                        watermarked_image.save(image_binary, format='PNG')
                        image_binary.seek(0)
                        await interaction.followup.send(file=discord.File(fp=image_binary, filename='watermarked_image.png'))
            except:
                await interaction.followup.send(traceback.format_exc())
        except:
            if arg is None:
                return
            await interaction.followup.send(f'{arg} is an unknown molecule.')

    @app_commands.command(name='smiles')
    async def smiles(self, interaction: discord.Interaction, molecules: str):
        try:
            args = shlex.split(molecules)
            try:
                for arg in args:
                    compounds = pcp.get_compounds(arg, 'name')
                    compound = compounds[0]
                    isomeric_smiles = compound.isomeric_smiles
                    await interaction.response.send_message(f'The isomeric SMILES for {arg} is: {isomeric_smiles}')
            except:
                await interaction.response.send_message(traceback.format_exc())
        except:
            if arg is None:
                return
            await interaction.response.send_message(f'{arg} is an unknown molecule.')

async def setup(bot: commands.Bot):
    await bot.add_cog(UserCog(bot))
