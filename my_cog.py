""" mycog.py
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

from discord.ext import commands
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import SimilarityMaps
from rdkit.DataStructs import FingerprintSimilarity, TanimotoSimilarity
from typing import Literal, Optional

import discord
import emoji as emoji_lib
import json
import io
import os
import pubchempy as pcp
import random
import unicodedata

class MyCog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.user_command_messages = {}

    @commands.Cog.listener()
    async def on_message(self, message: discord.Message) -> None:
        if message.author.bot:
            return
        # Record each command message
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
                mol = io.BytesIO()
                img.save(mol, format='PNG')
                mol.seek(0)
                file = discord.File(fp=mol, filename=f'Molecule.png')
                await ctx.send(file=file)
            else:
                mol = Chem.MolFromSmiles(args[0])
                if mol:
                     img = Draw.MolToImage(mol, size=(512, 512))
                     mol = io.BytesIO()
                     img.save(mol, format='PNG')
                     mol.seek(0)
                     file = discord.File(fp=mol, filename=f'Molecule.png')
                     await ctx.send(file=file)
                else:
                     await ctx.send('Invalid string. Enter a molecule name or SMILES string.')
                     await ctx.send(compounds)
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
                 await ctx.send(compounds)
                 await ctx.send('Invalid first string. Enter a molecule name or SMILES string.')
            if other_compounds:
                 other_compound_data = other_compounds[0].to_dict(properties=['isomeric_smiles'])
                 label_2 = other_compounds[0].synonyms[0]
                 refmol = Chem.MolFromSmiles(other_compound_data['isomeric_smiles'])
            else:
                 label_2 = 'refmol'
                 refmol = Chem.MolFromSmiles(args[1])
            if not refmol:
                await ctx.send(other_compounds)
                await ctx.send('Invalid second string. Enter a molecule name or SMILES string.')
            d2d = Draw.MolDraw2DCairo(512, 512)
            d2d2 = Draw.MolDraw2DCairo(512, 512)
            fp1 = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=256)
            fp2 = AllChem.GetMorganFingerprintAsBitVect(refmol, 2, nBits=256)
            similarity = TanimotoSimilarity(fp1, fp2)
            fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(refmol, mol, lambda m, i: SimilarityMaps.GetMorganFingerprint(m, i, radius=2, fpType='bv', nBits=2048), draw2d=d2d)
            d2d.FinishDrawing()
            buf = io.BytesIO(d2d.GetDrawingText())
            image1 = Image.open(buf)
            image1 = image1.convert("RGBA")
            fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(mol, refmol, lambda m, i: SimilarityMaps.GetMorganFingerprint(m, i, radius=2, fpType='bv', nBits=2048), draw2d=d2d2)
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
     #       for x in range(combined.width):
      #          for y in range(combined.height):
       #             pixel = combined.getpixel((x, y))
        #            if sum(pixel) >= 1020:
         #               draw.point((x, y), fill=(0, 0, 0, 0))
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
    async def emoji(self, ctx: commands.Context, emoji_character: str = None):
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

    @commands.command()
    async def how(self, ctx: commands.Context):
        help_message = """
**Bot Commands Help**

**!colorize <R> <G> <B>**
- **Description:** Changes your role color using RGB values.
- **Example:** `!colorize 255 0 0` (for red)

**!dmpurge**
- **Description:** Purges DM messages with Lucy.
- **Example:** `!dmpurge`

**!draw <molecule_name|SMILES> [<reference_molecule_name|SMILES>]**
- **Description:** Draws a molecule when provided a molecule name or SMILES, or compares two molecules graphically using RDKit.
- **Example:** `!draw C1CCCCC1` or `!draw LSD C1=CC=CC=C1`

**!emoji <emoji>**
- **Description:** Provides information about an emoji.
- **Example:** `!emoji 😀`

**!load <extension>**
- **Description:** Loads a cog.
- **Example:** `!load mycog`

**!purge**
- **Description:** Deletes up to 100 non-pinned messages from the channel.
- **Example:** `!purge`

**!reload <extension>**
- **Description:** Reloads a cog.
- **Example:** `!reload mycog`

**!smiles <chemical>**
- **Description:** Provides a SMILES code for a chemical.
- **Example:** `!smiles water`

**!unload <extension>**
- **Description:** Unloads a cog.
- **Example:** `!unload mycog`

**!sync [~|*|^] [<guild_id> ...]**
- **Description:** Synchronizes the bot's command tree.
  - `~`: Syncs commands to the current guild.
  - `*`: Copies global commands to the current guild and syncs.
  - `^`: Clears commands in the current guild.
  - Without arguments: Syncs commands globally.
  - With `guild_id` arguments: Syncs commands to the specified guilds.
- **Example:** `!sync` or `!sync ~` or `!sync *` or `!sync ^`

For additional assistance or if you encounter issues, please contact support.
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

    @commands.command(name='msds', help='Displays a shortened MSDS of a molecule.')
    async def msds(self, ctx, *, molecule_name: str):
        # Fetch molecule data from PubChem
        compounds = pcp.get_compounds(molecule_name, 'name')
        if not compounds:
            await ctx.send('Molecule not found.')
            return

        compound = compounds[0]
        compound_info = compound.to_dict(properties=['iupac_name', 'molecular_formula', 'molecular_weight'])
        # Fetch MSDS data
        msds_url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{compound.cid}"
        
        # Prepare MSDS message
        msds_message = (
            f"**Molecule Name:** {compound_info.get('iupac_name', 'N/A')}\n"
            f"**Molecular Formula:** {compound_info.get('molecular_formula', 'N/A')}\n"
            f"**Molecular Weight:** {compound_info.get('molecular_weight', 'N/A')} g/mol\n"
            f"**Description:** {compound_info.get('description', 'N/A')}\n"
            f"**MSDS URL:** [View MSDS]({msds_url})"
        )

        # Send the MSDS message
        await ctx.send(msds_message)

async def setup(bot: commands.Bot):
    await bot.add_cog(MyCog(bot))
