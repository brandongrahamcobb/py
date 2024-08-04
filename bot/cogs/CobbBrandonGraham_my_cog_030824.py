""" my_cog.py
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

from bs4 import BeautifulSoup
from bot.utils.helpers import load_config
from bot.utils.helpers import fetch_and_parse
from bot.utils.helpers import interact_with_chatgpt
from discord.ext import commands
from googleapiclient.discovery import build
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import SimilarityMaps
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
from rdkit.DataStructs import FingerprintSimilarity, TanimotoSimilarity
from typing import Literal, Optional

from rdkit.Chem import ChemicalFeatures
from rdkit import rdBase
from rdkit.RDPaths import RDDocsDir
from rdkit.RDPaths import RDDataDir
from rdkit.Chem.Draw import IPythonConsole

import discord
import emoji as emoji_lib
import json
import io
import os
import matplotlib.pyplot as plt
#import numpy
import pubchempy as pcp
#import random
import requests
#import traceback
import unicodedata


class MyCog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        config = load_config()
        self.user_command_messages = {}
        self.google_json_key = config.get('google_json_key')
        self.my_cse_id = '4231a573aa34241e4'
        self.headers = config.get('headers')

    def google_search(search_term, api_key, cse_id, **kwargs):
        service = build('customsearch', 'v1', developerKey=api_key)
        res = service.cse().list(q=search_term, cx=cse_id, **kwargs).execute()
        return res['items']

    @commands.Cog.listener()
    async def on_message(self, message: discord.Message) -> None:
        if message.author.bot:
            return
        if message.author != self.bot.user and message.content.startswith('!'):
            if message.author.id not in self.user_command_messages:
                self.user_command_messages[message.author.id] = []
            self.user_command_messages[message.author.id].append(message.id)

    @commands.Cog.listener()
    async def on_message_edit(self, before, after):
        if before.author == self.bot.user:
            return
        if after.content and after.content != before.content:
            if after.content.startswith('!'):
                ctx = await self.bot.get_context(after)
                if ctx.valid:
                    await self.bot.invoke(ctx)

    @commands.Cog.listener()
    async def on_error(self, event, *args, **kwargs) -> None:
        if event:
            logger.error(f'Error occurred in event {event}', exc_info=True)

    @commands.command(name='catch_all')
    async def catch_all(self, ctx, *, command_name: str):
        try:
            text = fetch_and_parse(f'https://www.cannabinoidswithturkey.com/cannabinoidinfosheets/')
            reftext = fetch_and_parse(f'https://www.cannabinoidswithturkey.com/cannabinoidinfosheets/{command_name}')
            if text.strip() != reftext.strip():
                await ctx.send(f'https://www.cannabinoidswithturkey.com/cannabinoidinfosheets/{command_name}')
            text = fetch_and_parse(f'https://www.cannabinoidswithturkey.com/info-sheets-na/')
            reftext = fetch_and_parse(f'https://www.cannabinoidswithturkey.com/info-sheets-na/{command_name}')
            if text.strip() != reftext.strip():
                await ctx.send(f'https://www.cannabinoidswithturkey.com/info-sheets-na/{command_name}')
            text = fetch_and_parse(f'https://legaldrugswithturkey.wixsite.com/legaldrugswithturkey/infosheets/')
            reftext = fetch_and_parse(f'https://legaldrugswithturkey.wixsite.com/legaldrugswithturkey/infosheets/{command_name}')
            if text.strip() != reftext.strip():
                await ctx.send(f'https://legaldrugswithturkey.wixsite.com/legaldrugswithturkey/infosheets/{command_name}')
        except Exception as e:
            await ctx.send(e)

    @commands.Cog.listener()
    async def on_command_error(self, ctx, error):
        if isinstance(error, commands.CommandNotFound):
            command_name = ctx.message.content[len('!'):].strip()
            catch_all_command = self.bot.get_command('catch_all')
            if catch_all_command:
                await ctx.invoke(catch_all_command, command_name=command_name)
        else:
            raise error

    @commands.command(description='')
    @commands.is_owner()
    async def colorize(self, ctx: commands.Context, *args):
        r = int(args[0])
        g = int(args[1])
        b = int(args[2])
        guildroles = await ctx.guild.fetch_roles()
        position = len(guildroles) - 1
        for arg in ctx.author.roles:
            if arg.name.isnumeric():
                await ctx.author.remove_roles(arg)
        for arg in guildroles:
            if arg.name.lower() == f'{r}{g}{b}':
                await ctx.author.add_roles(arg)
                await arg.edit(position=position)
                await ctx.send(f'I successfully changed your role color to {r}, {g}, {b}')
                return
        newrole = await ctx.guild.create_role(name=f'{r}{g}{b}', color=discord.Color.from_rgb(r, g, b), reason='new color')
        await newrole.edit(position=position)
        await ctx.author.add_roles(newrole)
        await ctx.send(f'I successfully changed your role color to {r}, {g}, {b}')

    @commands.command()
    async def dmpurge(self, ctx: commands.Context):
        async for message in ctx.author.history(limit=9999999):
            if message.author.id == self.bot.user.id:
                await message.delete()

    @commands.command(description='Compare molecules, fetch molecules.')
    async def draw(self, ctx: commands.Context, *args) -> None:
        compounds = pcp.get_compounds(args[0], 'name')
        if len(args) == 1:
            if compounds:
                compound_data = compounds[0].to_dict(properties=['isomeric_smiles'])
                mol = Chem.MolFromSmiles(compound_data['isomeric_smiles'])
                img = Draw.MolToImage(mol, size=(512, 512))
                img_bytes = io.BytesIO()
                img.save(img_bytes, format='PNG')
                img_bytes.seek(0)
                file = discord.File(fp=img_bytes, filename=f'Molecule.png')
                await ctx.send(file=file)
            else:
                mol = Chem.MolFromSmiles(args[0])
                if mol:
                     img = Draw.MolToImage(mol, size=(512, 512))
                     img_bytes = io.BytesIO()
                     img.save(img_bytes, format='PNG')
                     img_bytes.seek(0)
                     file = discord.File(fp=img_bytes, filename=f'Molecule.png')
                     await ctx.send(file=file)
                else:
                     embed = discord.Embed()
                     results = MyCog.google_search(args[0], self.google_json_key, self.my_cse_id, num=10)
                     for result in results:
                         embed.add_field(name=result['title'], value=result['link'], inline=False)
                     await ctx.send(embed=embed)
        elif len(args) == 2:
            other_compounds = pcp.get_compounds(args[1], 'name')
            if compounds:
                 compound_data = compounds[0].to_dict(properties=['isomeric_smiles'])
                 label_1 = compounds[0].synonyms[0]
                 mol = Chem.MolFromSmiles(compound_data['isomeric_smiles'])
            else:
                 label_1 = 'mol'
                 mol = Chem.MolFromSmiles(args[0])
            if not mol:
                 await ctx.send('Invalid first string. Enter a molecule name or SMILES string.')
            if other_compounds:
                 other_compound_data = other_compounds[0].to_dict(properties=['isomeric_smiles'])
                 label_2 = other_compounds[0].synonyms[0]
                 refmol = Chem.MolFromSmiles(other_compound_data['isomeric_smiles'])
            else:
                 label_2 = 'refmol'
                 refmol = Chem.MolFromSmiles(args[1])
            if not refmol:
                await ctx.send('Invalid second string. Enter a molecule name or SMILES string.')
            d2d = Draw.MolDraw2DCairo(512, 512)
            d2d2 = Draw.MolDraw2DCairo(512, 512)
            fp1 = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            fp2 = AllChem.GetMorganFingerprintAsBitVect(refmol, 2, nBits=2048)
            similarity = TanimotoSimilarity(fp1, fp2)
            fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(refmol, mol, lambda m, i: SimilarityMaps.GetMorganFingerprint(m, i, radius=2, fpType='bv', nBits=8192), draw2d=d2d)
            d2d.FinishDrawing()
            buf = io.BytesIO(d2d.GetDrawingText())
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
 #           for x in range(combined.width):
  #              for y in range(combined.height):
   #                 pixel = combined.getpixel((x, y))
    #                if sum(pixel) >= 1020:
     #                   draw.point((x, y), fill=(0, 0, 0, 0))
            font = ImageFont.load_default(40)
            draw.text(((width - draw.textlength(label_1, font=font)) / 2, 10), label_1, fill="black", font=font)
            draw.text((width + (width - draw.textlength(label_2, font=font)) / 2, 10), label_2, fill="black", font=font)
            text = f"Tanimoto Similarity: {similarity:.5f}"
            draw.text(((combined.width - draw.textlength(text, font=font)) / 2, height - 60), text, fill="black", font=font)
            output_buffer = io.BytesIO()
            combined.save(output_buffer, format='PNG')
            output_buffer.seek(0)
            file = discord.File(fp=output_buffer, filename=f'Molecule.png')
            await ctx.send(file=file)

    @commands.command()
    async def how(self, ctx: commands.Context):
        help_message = """
        """
        await ctx.send(help_message)

    @commands.command(name='load', hidden=True)
    @commands.has_permissions(ban_members=True)
    async def load(self, ctx: commands.Context, *, module: str):
        try:
            await ctx.bot.load_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')

    @commands.command()
    async def lucy(self, ctx):
         file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'my_cog.py')
         file = discord.File(file_path, filename='my_cog.py')
         await ctx.send("Here is your file:", file=file)

    @commands.command()
    @commands.has_permissions(ban_members=True)
    async def reload(self, ctx: commands.Context, *, module: str):
        try:
            await ctx.bot.reload_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')

    @commands.command()
    async def smiles(self, ctx: commands.Context, *, chemical_name: str):
        compounds = pcp.get_compounds(chemical_name, 'name')
        if not compounds:
            await ctx.send(f"No compound found for the name '{chemical_name}'")
            return
        compound = compounds[0]
        canonical_smiles = compound.canonical_smiles
        await ctx.send(f"The canonical SMILES for '{chemical_name}' is: {canonical_smiles}")

    @commands.command(name='sync', hidden=True)
    @commands.has_permissions(ban_members=True)
    async def sync(self, ctx: commands.Context, guilds: commands.Greedy[discord.Object], spec: Optional[Literal["~", "*", "^"]] = None) -> None:
        if not guilds:
            if spec == "~":
                synced = await ctx.bot.tree.sync(guild=ctx.guild)
            elif spec == "*":
                ctx.bot.tree.copy_global_to(guild=ctx.guild)
                synced = await ctx.bot.tree.sync(guild=ctx.guild)
            elif spec == "^":
                ctx.bot.tree.clear_commands(guild=ctx.guild)
                await ctx.bot.tree.sync(guild=ctx.guild)
                synced = []
            else:
                synced = await ctx.bot.tree.sync()
            await ctx.send(
                f"Synced {len(synced)} commands {'globally' if spec is None else 'to the current guild.'}"
            )
            return
        ret = 0
        for guild in guilds:
            try:
                await ctx.bot.tree.sync(guild=guild)
            except discord.HTTPException:
                pass
            else:
                ret += 1
        await ctx.send(f"Synced the tree to {ret}/{len(guilds)}.")

    @commands.command(name='unload', hidden=True)
    @commands.has_permissions(ban_members=True)
    async def unload(self, ctx: commands.Context, *, module: str):
        try:
            await ctx.bot.unload_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')

    @commands.command()
    @commands.has_permissions(manage_messages=True)
    async def purge(self, ctx, option=None, limit: int = 100):
        # Check if the limit is valid
        if limit <= 0 or limit > 100:
            await ctx.send("Please provide a limit between 1 and 100.")
            return
        if option == "bot":
            # Purge only bot messages
            def is_bot_message(message):
                return message.author == self.bot.user
            deleted = await self.purge_messages(ctx, limit, is_bot_message)
            await ctx.send(f"Deleted {deleted} bot messages.")
        elif option == "all":
            # Purge all messages
            deleted = await self.purge_messages(ctx, limit)
            await ctx.send(f"Deleted {deleted} messages.")
        elif option == "user":
            # Purge user's messages
            def is_user_message(message):
                return message.author == ctx.author
            deleted = await self.purge_messages(ctx, limit, is_user_message)
            await ctx.send(f"Deleted {deleted} of your messages.")
        elif option == "commands":
            # Delete all commands from the user
            if ctx.author.id in self.user_command_messages:
                message_ids = self.user_command_messages[ctx.author.id]
                deleted = 0
                for message_id in message_ids:
                    try:
                        message = await ctx.channel.fetch_message(message_id)
                        await message.delete()
                        deleted += 1
                    except discord.NotFound:
                        continue
                del self.user_command_messages[ctx.author.id]
                await ctx.send(f"Deleted {deleted} of your commands.")
            else:
                await ctx.send("No commands found to delete.")
        else:
            await ctx.send("Invalid option. Use `bot`, `all`, `user`, or `commands`.")

    async def purge_messages(self, ctx, limit, check=None):
        """Purge messages with an optional check."""
        deleted = 0
        async for message in ctx.channel.history(limit=limit):
            if check is None or check(message):
                await message.delete()
                deleted += 1
        return deleted

    @commands.command(name="info")
    async def info(self, ctx, info_type: str = None, *, argument: str = None):
        if info_type is None:
            await ctx.send("Please provide the type of information you need. Options are: `emoji`, `msds`.")
            return
        if info_type.lower() == "emoji":
            if argument is None:
                await ctx.send("Please provide the Unicode emoji character.")
            else:
                await self.get_emoji_info(ctx, argument)
        elif info_type.lower() == "msds":
            if argument is None:
                await ctx.send("Please provide the molecule name.")
            else:
                await self.get_msds_info(ctx, argument)
        else:
            await ctx.send("Invalid info type. Options are: `emoji`, `msds`.")

    async def get_emoji_info(self, ctx, emoji_character):
        if emoji_character is None:
            await ctx.send('Please provide a Unicode emoji character.')
            return
        unicode_name = emoji_lib.demojize(emoji_character)
        if unicode_name.startswith(':') and unicode_name.endswith(':'):
            unicode_name = unicode_name[1:-1]
            code_points = ' '.join(f'U+{ord(c):04X}' for c in emoji_character)
            description = unicodedata.name(emoji_character, 'No description available')
            unicode_info = (
                f'**Unicode Emoji Character:** {emoji_character}\n'
                f'**Unicode Emoji Name:** {unicode_name}\n'
                f'**Code Points:** {code_points}\n'
                f'**Description:** {description}'
            )
            await ctx.send(unicode_info)
        else:
            await ctx.send('Unicode emoji not found. Make sure it is a valid Unicode emoji character.')

    async def get_msds_info(self, ctx, molecule_name: str):
        # Fetch molecule data from PubChem
        compounds = pcp.get_compounds(molecule_name, 'name')
        if not compounds:
            await ctx.send('Molecule not found.')
            return
        compound = compounds[0]
        compound_info = compound.to_dict(properties=['iupac_name', 'molecular_formula', 'molecular_weight', 'synonyms'])
        description = compound_info.get('synonyms', ['No description available'])[0]
        msds_url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{compound.cid}"
        fda_guide_url = self.get_fda_prescriber_guide(molecule_name)
        msds_message = (
            f"**Molecule Name:** {compound_info.get('iupac_name', 'N/A')}\n"
            f"**Molecular Formula:** {compound_info.get('molecular_formula', 'N/A')}\n"
            f"**Molecular Weight:** {compound_info.get('molecular_weight', 'N/A')} g/mol\n"
            f"**Description:** {description}\n"
            f"**MSDS URL:** [View MSDS]({msds_url})\n"
            f"**FDA Prescriber Guide:** {fda_guide_url if fda_guide_url else 'No prescriber guide found.'}"
        )
        await ctx.send(msds_message)

    def get_fda_prescriber_guide(self, molecule_name: str) -> str:
        search_query = f"{molecule_name} prescriber guide"
        google_search_url = f"https://www.google.com/search?q={search_query}"
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
        }
        response = requests.get(google_search_url, headers=headers)
        if response.status_code != 200:
            return None
        soup = BeautifulSoup(response.text, "html.parser")
        for result in soup.find_all("a"):
            href = result.get("href")
            if href and "fda.gov" in href:
                return href.split("&")[0].replace("/url?q=", "")
        return None

    @commands.command(name='botinfo')
    @commands.is_owner()  # Ensures only the bot owner can use this command
    async def botinfo(self, ctx, member: discord.Member, prompt: str):
        """DMs the specified user with the bot's information."""
        bot_info = (
            "Here is some information about the bot:\n"
            f"**Name:** {self.bot.user.name}\n"
            f"**ID:** {self.bot.user.id}\n"
            f"**Description:** This is a bot created by Brandon Graham Cobb for various functionalities.\n"
            f"**Invite Link:** [Invite Bot](https://discord.com/oauth2/authorize?client_id=302202228016414721)\n"
            f"**Support Server:** [Support Server](https://discord.gg/dNfUn8MeYN)\n"
            f"**P.S.:** {interact_with_chatgpt(prompt)}"
        )
        try:
            await member.send(bot_info)
            await ctx.send(f"Sent bot information to {member.mention}.")
        except discord.Forbidden:
            await ctx.send(f"Failed to send DM to {member.mention}. They might have DMs disabled.")

async def setup(bot: commands.Bot):
    await bot.add_cog(MyCog(bot))
