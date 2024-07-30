
# cogs/vegan_university.py

from bot import Lucy
from datetime import datetime
from discord.ext import commands
from member_data import MemberData
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import SimilarityMaps
from rdkit.DataStructs import FingerprintSimilarity, TanimotoSimilarity
from typing import Literal, Optional

import discord
import emoji as emoji_lib
import io
import os
import pubchempy as pcp
import unicodedata

class VeganUniversity(commands.Cog):

    def __init__(self, bot: Lucy):
        self.bot: Lucy = bot
        self.bot.member_data = {}

    def version():
        if not os.path.exists('../txt/version.txt'):
             with open('../txt/version.txt', 'w') as f:
                 f.write('1.0.0')
        with open('../txt/version.txt', 'r') as f:
             version = f.read().strip()
        major, minor, patch = map(int, version.split('.'))
        patch += 1
        if patch >= 10:
            patch = 0
            minor += 1
        if minor >= 10:
            minor = 0
            major += 1
        version = f"{major}.{minor}.{patch}"
        with open('../txt/version.txt', 'w') as f:
            f.write(version)
        return version

    @commands.Cog.listener()
    async def on_ready(self):
        version = VeganUniversity.version()
        print('Logged on as', self.bot.user)
        print(f'Bot ID: {self.bot.user.id}')
        print(f'Bot Version: {version}')
        print(f'Connected to {len(self.bot.guilds)} guilds:')
        for guild in self.bot.guilds:
            print(f' - {guild.name} (id: {guild.id})')

    @commands.Cog.listener()
    async def on_member_join(self, member: discord.Member):
        self.bot.member_data[member.id] = MemberData(join_time=datetime.utcnow())
        print(f'{member} has joined the server.')

    @commands.Cog.listener()
    async def on_member_remove(self, member: discord.Member):
        if member.id in self.bot.member_data:
            del self.bot.member_data[member.id]
        print(f'{member} has left the server.')

    @commands.Cog.listener()
    async def on_command(self, ctx: commands.Context):
        if ctx.author.id in self.bot.member_data:
            exp_points = len(ctx.message.content) - 1  # Calculate experience based on command length
            self.member_data[ctx.author.id].increment_command(ctx.command.name, exp_points)
            print(f'{ctx.author} used {ctx.command.name} and gained {exp_points} XP.')

    @commands.Cog.listener()
    async def on_message(self, message: discord.Message) -> None:
        if message.author.bot:
            return

    @commands.Cog.listener()
    async def on_error(event, *args, **kwargs) -> None:
        if event:
             return

    @commands.command()
    @commands.is_owner()
    async def purge(self, ctx: commands.Context):
        def not_pinned(msg):
            return not msg.pinned
        purged = await ctx.channel.purge(limit=100, check=not_pinned)
        await ctx.send(f"Successfully removed {len(purged)} non-pinned messages!")

    @commands.command()
    @commands.is_owner()
    async def dmpurge(self, ctx: commands.Context):
        async for message in self.bot.get_user(154749533429956608).history(limit=9999999):
            if message.author.id == self.bot.user.id:
                await message.delete()

    @commands.command(name='load', hidden=True)
    @commands.has_permissions(ban_members=True)
    async def load(self, ctx: commands.Context, *, module: str):
        try:
            await ctx.bot.load_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')

    @commands.command(name='reload', hidden=True)
    @commands.has_permissions(ban_members=True)
    async def reload(self, ctx: commands.Context, *, module: str):
        try:
            await ctx.bot.reload_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')

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
        for x in range(combined.width):
            for y in range(combined.height):
                pixel = combined.getpixel((x, y))
                if sum(pixel) >= 1020:
                    draw.point((x, y), fill=(0, 0, 0, 0))
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
        help_message = (
            "**Bot Commands**\n\n"
            "**!colorize**: Change your nickname color. Only the bot owner can use this command.\n"
            "**Usage**: `!colorize <R> <G> <B>`\n\n"
            "**!draw**: Draws a molecule by SMILES or compares two molecules displaying a Tanimoto similarity along with their 2D images.\n"
            "**Usage**: `!draw <SMILES> or !draw <MOLECULE> <REFERENCEMOLECULE>`\n\n"
            "**!emoji**: Provides information about a given Unicode emoji character, including its Unicode name, code points, and description.\n"
            "**Usage**: `!emoji <emoji_character>`\n\n"
            "**!penji**: Connects you to the MCC tutoring service.\n"
            "**Usage**: `!penji`\n\n"
            "**!power**: Compares two molecules based on their SMILES strings and displays their Tanimoto similarity along with their 2D images.\n"
            "**Usage**: `!power <SMILES1> <TITLE1> <SMILES2> <TITLE2>`\n\n"
            "**!purge**: Deletes up to 100 non-pinned messages from the channel. Only the bot owner can use this command.\n"
            "**Usage**: `!purge`\n\n"
            "**!smiles**: Retrieves the SMILES of any molecule in the PubChem database..\n"
            "**Usage**: `!smiles <MOLECULE>`"
        )
        await ctx.send(help_message)

    @commands.command()
    async def level(self, ctx: commands.Context):
        if ctx.author.id in self.bot.member_data:
            member_data = self.bot.member_data[ctx.author.id]
            level = member_data.get_current_level()
            await ctx.send(f'You are currently level {level}.')
        else:
            await ctx.send('Member data not found.')

    @commands.command()
    async def penji(self, ctx: commands.Context, *args) -> None:
        await ctx.send('Sign into https://cutt.ly/mccpenjikiosk-cc')

    @commands.command()
    async def power(self, ctx: commands.Context, *args):
        mol = Chem.MolFromSmiles(args[0])
        refmol = Chem.MolFromSmiles(args[2])
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
        for x in range(combined.width):
            for y in range(combined.height):
                pixel = combined.getpixel((x, y))
                if sum(pixel) >= 1020:
                    draw.point((x, y), fill=(0, 0, 0, 0))
        font = ImageFont.load_default(40)
        text = args[1]
        draw.text(((width - draw.textlength(text, font=font)) / 2, 10), text, fill="black", font=font)
        text = args[3]
        draw.text((width + (width - draw.textlength(text, font=font)) / 2, 10), text, fill="black", font=font)
        text = f"Tanimoto Similarity: {similarity:.5f}"
        draw.text(((combined.width - draw.textlength(text, font=font)) / 2, height - 60), text, fill="black", font=font)
        output_buffer = io.BytesIO()
        combined.save(output_buffer, format='PNG')
        output_buffer.seek(0)
        file = discord.File(fp=output_buffer, filename=f'Molecule.png')
        await ctx.send(file=file)

    @commands.command()
    async def smiles(self, ctx: commands.Context, *, chemical_name: str):
        compounds = pcp.get_compounds(chemical_name, 'name')
        if not compounds:
            await ctx.send(f"No compound found for the name '{chemical_name}'")
            return
        compound = compounds[0]
        canonical_smiles = compound.canonical_smiles
        await ctx.send(f"The canonical SMILES for '{chemical_name}' is: {canonical_smiles}")

    @commands.command()
    async def version(self, ctx: commands.Context):
        version = Lucy.version()
        await ctx.send(f'Bot Version: {version}')

async def setup(bot: Lucy):
    await bot.add_cog(VeganUniversity(bot))
