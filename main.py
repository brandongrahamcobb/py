from datetime import datetime
from discord.ext import commands
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import SimilarityMaps
from rdkit.DataStructs import FingerprintSimilarity, TanimotoSimilarity
from typing import Literal, Optional

import asyncio
import discord
import emoji as emoji_lib
import json
import io
import logging
import logging.handlers
import os
import pubchempy as pcp
import random
import requests
import subprocess
import sys
import threading
import unicodedata

xp_table = [1, 15, 34, 57, 92, 135, 372, 560, 840, 1242, 1144, 1573, 2144, 2800, 3640, 4700, 5893, 7360, 9144, 11120, 
          13477, 16268, 19320, 22880, 27008, 31477, 36600, 42444, 48720, 55813, 63800, 86784, 98208, 110932, 124432, 
          139372, 155865, 173280, 192400, 213345, 235372, 259392, 285532, 312928, 342624, 374760, 408336, 445544, 
          483532, 524160, 567772, 598886, 631704, 666321, 702836, 741351, 781976, 824828, 870028, 917625, 967995, 
          1021041, 1076994, 1136013, 1198266, 1263930, 1333194, 1406252, 1483314, 1564600, 1650340, 1740778, 1836173, 
          1936794, 2042930, 2154882, 2272970, 2397528, 2528912, 2667496, 2813674, 2967863, 3130502, 3302053, 3483005, 
          3673873, 3875201, 4087562, 4311559, 4547832, 4797053, 5059931, 5337215, 5629694, 5938202, 6263614, 6606860, 
          6968915, 7350811, 7753635, 8178534, 8626718, 9099462, 9598112, 10124088, 10678888, 11264090, 11881362, 
          12532461, 13219239, 13943653, 14707765, 15513750, 16363902, 17260644, 18206527, 19204245, 20256637, 21366700, 
          22537594, 23772654, 25075395, 26449526, 27898960, 29427822, 31040466, 32741483, 34535716, 36428273, 38424542, 
          40530206, 42751262, 45094030, 47565183, 50171755, 52921167, 55821246, 58880250, 62106888, 65510344, 69100311, 
          72887008, 76881216, 81094306, 85594273, 90225770, 95170142, 100385466, 105886589, 111689174, 117809740, 
          124265714, 131075474, 138258410, 145834970, 153826726, 162256430, 171148082, 180526997, 190419876, 200854885, 
          211861732, 223471711, 223471711, 248635353, 262260570, 276632449, 291791906, 307782102, 324648562, 342439302, 
          361204976, 380999008, 401877754, 423900654, 447130410, 471633156, 497478653, 524740482, 553496261, 583827855, 
          615821622, 649568646, 685165008, 722712050, 762316670, 804091623, 848155844, 894634784, 943660770, 995373379, 
          1049919840, 1107455447, 1168144006, 1232158297, 1299680571, 1370903066, 1446028554, 1525246918, 1608855764, 
          1697021059, 2000000000]

# Utility functions
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

def verify_discord_token(token):
    headers = {
        "Authorization": f"Bot {token}"
    }
    response = requests.get("https://discord.com/api/v10/users/@me", headers=headers)
    return response.status_code == 200

def get_level_from_xp(xp):
    """ Calculate the level based on XP. """
    level = 0
    while level < len(xp_table) and xp >= xp_table[level]:
        level += 1
    return level

def is_root():
    return os.geteuid() == 0  # Check if running as root

def load_data():
    try:
        with open('xp_data.json', 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        return {}

def save_data(data):
    with open('xp_data.json', 'w') as f:
        json.dump(data, f, indent=4)

# Setup logging
logger = logging.getLogger('discord')
logger.setLevel(logging.INFO)
handler = logging.handlers.RotatingFileHandler(
    filename=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'log', 'discord.log'),
    encoding='utf-8',
    maxBytes=32 * 1024 * 1024,
    backupCount=5,
)
dt_fmt = '%Y-%m-%d %H:%M:%S'
formatter = logging.Formatter('[{asctime}] [{levelname:<8} {name}: {message}', dt_fmt, style='{')
handler.setFormatter(formatter)
logger.addHandler(handler)

# Bot setup
intents = discord.Intents.all()
intents.message_content = True
exts = []  # Add any cogs here if needed

user_data = load_data()

class Lucy(commands.Bot):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    async def on_ready(self):
        logger.info(f'Logged in as {self.user.name} ({self.user.id})')
        logger.info(f'Bot version: {version()}')
        logger.info(f'Connected to {len(self.guilds)} guilds:')
        for guild in self.guilds:
            logger.info(f' - {guild.name} (id: {guild.id})')

class ScriptCog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.member_data = {}

    @commands.Cog.listener()
    async def on_member_join(member):
        user_id = str(member.id)
        if user_id not in user_data:
            user_data[user_id] = {'xp': 0, 'level': 0}
            save_data(user_data)
            print(f"Initialized XP data for new member: {member.name}")

    @commands.Cog.listener()
    async def on_member_remove(member):
        user_id = str(member.id)
        if user_id in user_data:
            del user_data[user_id]
            save_data(user_data)
            print(f"Reset XP data for member: {member.name}")

    @commands.Cog.listener()
    async def on_message(self, message: discord.Message) -> None:
        if message.author.bot:
            return
        user_id = str(message.author.id)
        if user_id not in user_data:
            user_data[user_id] = {'xp': 0, 'level': 0}
        current_xp = user_data[user_id]['xp']
        current_level = user_data[user_id]['level']
        if current_level < len(xp_table):
            max_increase = (xp_table[current_level] - current_xp) // 10
        else:
            max_increase = 50
        xp_increase = random.randint(1, max_increase)
        user_data[user_id]['xp'] += xp_increase
        new_level = get_level_from_xp(user_data[user_id]['xp'])
        user_data[user_id]['level'] = new_level
        save_data(user_data)
        await bot.process_commands(message)

    @commands.Cog.listener()
    async def on_error(self, event, *args, **kwargs) -> None:
        if event:
            logger.error(f'Error occurred in event {event}', exc_info=True)

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

    @commands.command()
    @commands.has_permissions(administrator=True)
    async def addxp(ctx, amount: int, member: discord.Member = None):
        member = member or ctx.author
        user_id = str(member.id)
        if user_id not in user_data:
             user_data[user_id] = {'xp': 0, 'level': 0}
        user_data[user_id]['xp'] += amount
        user_data[user_id]['level'] = get_level_from_xp(user_data[user_id]['xp'])
        await ctx.send(f'{amount} XP added to {member.mention}.')
        save_data(user_data)
        await xp(ctx)

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
    async def lucy(self, ctx):
         # Path to the file you want to send
         file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'main.py')
         # Create a discord.File object
         file = discord.File(file_path, filename='main.py')
         # Send the file to the channel
         await ctx.send("Here is your file:", file=file)

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
        version = version()
        await ctx.send(f'Bot version: {version}')

    @commands.command()
    async def xp(self, ctx):
        user_id = str(ctx.author.id)
        if user_id not in user_data:
            user_data[user_id] = {'xp': 0, 'level': 0}
        xp = user_data[user_id]['xp']
        level = user_data[user_id]['level']
        await ctx.send(f'{ctx.author.mention} - Level: {level}, XP: {xp}')
        user_id = str(ctx.author.id)
        if user_id not in user_data:
            user_data[user_id] = {'xp': 0, 'level': 0}

async def setup():
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Welcome, to initial setup~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    system = 'Linux' if os.name != 'nt' else 'Windows'
    
    if system == 'Linux':
        if is_root():
            print('Running as root, updating system...')
            subprocess.run(['pacman', '-Syu'], check=True)
        else:
            print('Not running as root, skipping system update.')

    dir = os.path.dirname(os.path.abspath(__file__))
    updir = os.path.join(dir, '..')
    dir_venv = os.path.join(updir, 'activate')

    # Create virtual environment
    subprocess.run([sys.executable, '-m', 'venv', dir_venv], check=True)
    
    requirements = ["asyncpraw", "discord.py", "emoji", "pubchempy", "rdkit", "pillow", "requests"]
    py = os.path.join(dir_venv, 'bin', 'python') if system == 'Linux' else sys.executable
    subprocess.run([py, '-m', 'pip', 'install'] + requirements, check=True)

    # Create necessary directories
    dir_json = os.path.join(updir, 'json')
    dir_log = os.path.join(updir, 'log')
    dir_txt = os.path.join(updir, 'txt')

    os.makedirs(dir_json, exist_ok=True)
    os.makedirs(dir_log, exist_ok=True)
    os.makedirs(dir_txt, exist_ok=True)

    print(f'System: {system}')
    print(f'Location: {dir}')
    print(f'Python Location: {sys.executable}')
    print(f'Python Version: {sys.version}')
    
    file_json = os.path.join(dir_json, 'config.json')

    # Check if config.json already exists
    if os.path.exists(file_json):
        with open(file_json, 'r') as f:
            config = json.load(f)
        if 'token' in config and verify_discord_token(config['token']):
            print('Token is valid and already configured.')
            return  # Skip token input if it's valid and exists

    # Prompt for token if not present or invalid
    while True:
        token = input('Enter your bot token: ')
        if verify_discord_token(token):
            print('Configuring....')
            with open(file_json, 'w') as f:
                json.dump({
                    'token': token,
                    'prefix': "!",
                    'os': system
                }, f, indent=4)
            break
        else:
            print('Enter a valid token.')

async def main():
    await setup()
    
    dir = os.path.dirname(os.path.abspath(__file__))
    updir = os.path.join(dir, '..')
    file_json = os.path.join(updir, 'json', 'config.json')

    with open(file_json, 'r') as f:
        data = json.load(f)
        token = data['token']

    bot = Lucy(
        command_prefix=commands.when_mentioned_or('!'),
        intents=intents,
    )
    await bot.add_cog(ScriptCog(bot))  # Add the cog to the bot

    try:
        await bot.start(token)
    except Exception as e:
        logger.error('Failed to start the bot.', exc_info=e)

# Run the setup and bot
if __name__ == '__main__':
    asyncio.run(main())
