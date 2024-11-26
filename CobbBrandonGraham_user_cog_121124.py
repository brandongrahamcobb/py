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

from bot.main import Spooky
from bs4 import BeautifulSoup
from discord import app_commands, Embed
from discord.ext import commands
from rdkit import Chem
from rdkit.Chem import Crippen
from typing import List

import bot.utils.helpers as helpers

import discord
import bs4
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
        self.config = bot.config
        self.stacks = {}
        self.diaminobenzidine = helpers.get_mol('3,3-diaminobenzidine')
        self.benadryl = helpers.get_mol('diphenhydramine')

    def get_user_stack(self, user_id: int):
        return self.stacks.get(user_id, [])

    def set_user_stack(self, user_id: int, molecule_names: List[str]):
        mol_objects = [helpers.get_mol(name) for name in molecule_names]
        self.stacks[user_id] = mol_objects

    @commands.command(name='stack', description='Chose your guide.')
    async def analog(self, ctx: commands.Context, *args):
        new_stack = list(args)
        if args:
            self.set_user_stack(user_id=ctx.author.id, molecule_names=args)
        while True:
            await ctx.send(f'{ctx.author.name}\'s Comparators:')
#            if self.get_user_stack(ctx.author.id):
            await ctx.send([helpers.get_molecule_name(mol) for mol in self.get_user_stack(ctx.author.id)])
            await ctx.send('Which molecule would you like to interview?')
            response = await self.bot.wait_for(
                'message',
                timeout=20000.0,
                check=lambda message: message.author == ctx.author and message.channel == ctx.channel
            )
#            original = max(self.get_user_stack(ctx.author.id), key=lambda mol: helpers.get_proximity(mol, self.nicotine))
            try:
                new = helpers.get_mol(response.content)
                molecules = []
                proximities = []
                for arg in args:
                    molecules.append(helpers.get_mol(arg))
                for arg in molecules:
                    proximities.append(helpers.get_proximity(new, arg))
                    
#                original_proximity = helpers.get_proximity(new, original)
 #               nicotine_proximity = helpers.get_proximity(new, self.nicotine)
#                diaminobenzidine_proximity = helpers.get_proximity(self.diaminobenzidine, new)
 #               benadryl_proximity = helpers.get_proximity(self.benadryl, new)
  #              if original_proximity > diaminobenzidine_proximity:
   #                 new_stack.append(response.content)
    #                embed = helpers.get_sds(response.content)
     #               await ctx.send(embed=embed)
      #              self.set_user_stack(user_id=ctx.author.id, molecule_names=new_stack)
#                else:
 #                   await ctx.send(f'DANGER: {response.content}')
                #await ctx.send(f'{(100 * nicotine_proximity):.3f}% {response.content} {(100 * diaminobenzidine_proximity):.3f}% {response.content} {(100 * lsd_proximity):.3f}%')
                #await ctx.send(f'Benadryl: {(100 * benadryl_proximity):.3f}%\n Diaminobenzidine: {(100 * diaminobenzidine_proximity):.3f}% You: {response.content}')
                embed = Embed(title='Tanimoto Table', color=0x00ff00)
                for molecule, proximity in zip(args, proximities):
                    embed.add_field(name='Molecule', value=molecule, inline=True)
                    embed.add_field(name="Similarity", value=f'{proximity:.3f}', inline=True)
                    embed.add_field(name="\u200b", value='\u200b', inline=True)
                string_builder = ''
                for field in embed.fields:
                     string_builder += f'{field.name}: {field.value}\n'
                await ctx.send(embed=embed)
                #await ctx.send(f'Benadryl: {(100 * benadryl_proximity):.3f}%\n Diaminobenzidine: {(100 * diaminobenzidine_proximity):.3f}% You: {response.content}')
            except Exception as e:
                await ctx.send(e)
                break

    @commands.hybrid_command(name='compare', description='Draw a molecule or compare molecules by their names.')
    async def compare(self, ctx: commands.Context, *, molecules: str) -> None:
        if ctx.interaction:
            await ctx.interaction.response.defer(ephemeral=True)
        args = shlex.split(molecules)
        mol1 = helpers.get_mol(args[0])
        mol2 = helpers.get_mol(args[1])
        image = helpers.view_difference(mol1, mol2)
        await ctx.send(file=image)

    @commands.hybrid_command(name='sim')
    async def sim(self, ctx: commands.Context, *, molecules: str):
        args = shlex.split(molecules)
        similarity = helpers.get_proximity(helpers.get_mol(args[0]), helpers.get_mol(args[1]))
        await ctx.send(similarity)

    @commands.hybrid_command(name='logp')
    async def logp(self, ctx: commands.Context, *, molecules: str):
        args = shlex.split(molecules)
        for arg in args:
            compounds = pcp.get_compounds(arg, 'name')
            compound_data = compounds[0].to_dict(properties=['isomeric_smiles'])
            mol = Chem.MolFromSmiles(compound_data['isomeric_smiles'])
            log_p = Crippen.MolLogP(mol)
            await ctx.send(f'Your octanol:water coefficient is: {log_p}')

    @commands.hybrid_command(name='msds')
    async def msds(self, ctx: commands.Context, *, argument: str):
       embed = helpers.get_sds(argument)
       await ctx.send(embed=embed)

    @commands.hybrid_command(name='smiles')
    async def smiles(self, ctx: commands.Context, *, molecules: str):
        try:
            if ctx.interaction:
                 await ctx.interaction.response.defer(ephemeral=True)
            args = shlex.split(molecules)
            try:
                for arg in args:
                    compounds = pcp.get_compounds(arg, 'name')
                    compound = compounds[0]
                    isomeric_smiles = compound.isomeric_smiles
                    await ctx.send(f'The isomeric SMILES for {arg} is: {isomeric_smiles}')
            except:
                await ctx.send(traceback.format_exc())
        except:
            if arg is None:
                return
            await ctx.send(f'{arg} is an unknown molecule.')

async def setup(bot: commands.Bot):
    await bot.add_cog(UserCog(bot))
