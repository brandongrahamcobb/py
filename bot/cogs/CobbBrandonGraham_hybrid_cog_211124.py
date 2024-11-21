''' hybrid.py The purpose of this program is to provide the core functionality to Vyrtuous.
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
import shlex
import traceback

class Hybrid(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.indica = self.bot.get_cog('Indica')
        self.sativa = self.bot.get_cog('Sativa')

    @commands.command(description='')
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

#    @commands.command(name='colors')
#    async def colors(self, ctx: commands.Context, *args):
#        try:
#            attachment = ctx.message.attachments[0]
#            img_data = requests.get(attachment.url).content
#            image = Image.open(io.BytesIO(img_data))
#            image = image.resize((100, 100))
#            image = image.convert('RGB')
#            pixels = list(image.getdata())
#            pixels = [pixel for pixel in pixels if not (pixel[0] > 150 and pixel[1] > 150 and pixel[2] > 150)]
#            pixels = [pixel for pixel in pixels if not (pixel[0] < 10 and pixel[1] < 10 and pixel[2] < 10)]
#            color_counts = Counter(pixels)
#            predominant_colors = color_counts.most_common(int(args[0]))
#            message = "Predominant colors (excluding whites above (150, 150, 150)):\n"
#            for color, count in predominant_colors:
#                message += f"Color: {color} Count: {count}\n"
#            await ctx.send(message)
#        except Exception as e:
#            await ctx.send(e)

    @commands.hybrid_command(name='get')
    async def get(self, ctx: commands.Context, *, argument: str):
        try:
            if ctx.interaction:
                await ctx.interaction.response.defer(ephemeral=True)
            file = helpers.stable_cascade(argument)
            if isinstance(file, discord.File):
                await ctx.send(file=file)
            else:
                await ctx.send(f"Error generating image: {file}")
        except Exception as e:
            await ctx.send(traceback.format_exc())

    @commands.command(name='load', hidden=True)
    async def load(self, ctx: commands.Context, *, module: str):
        try:
            await ctx.bot.load_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')

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
                await ctx.send('Invalid option. Use `compare`, `glow`, `gsrs`, or `shadow`.')
        except Exception as e:
            print(f'An error occurred: {traceback.format_exc()}')

    @commands.hybrid_command(name='recipe', description='Simplify a recipe.')
    async def recipe(self, ctx: commands.Context, url=None):
        try:
            if ctx.interaction:
                await ctx.interaction.response.defer(ephemeral=True)
            if not url:
                await ctx.send('No recipe url provided.')
            embed = helpers.get_recipe(url)
            await ctx.send(embed)
        except Exception as e:
            await ctx.send(e)

    @commands.hybrid_command()
    async def reload(self, ctx: commands.Context, *, module: str):
        try:
            await ctx.bot.reload_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')

    @commands.command(name='script', description='Usage !script <NIV/ESV> <Book>.<Chapter>.<Verse>')
    async def script(self, ctx: commands.Context, version: str, *, reference: str):
         await ctx.send(helpers.get_scripture(version, reference))

    @commands.hybrid_command(name='warn',description='removes another user\'s message')
    async def warn(self, ctx: commands.Context, ):
        try:
            message = await ctx.channel.fetch_message(message.id)
            if message:
                spoiler_content = f"||{message.content}||"
                await message.delete()
                await ctx.send(f"{message.author}: {spoiler_content}")
            else:
                await ctx.send("Message not found.")
        except Exception as e:
            await ctx.send(f'An error occurred: {e}')

    @commands.command(name='vegan')
    async def vegan(self, ctx: commands.Context, *args):
        try:
            await ctx.send('https://linktr.ee/veganresource')
        except:
            await ctx.send('Fail')

async def setup(bot: commands.Bot):
    await bot.add_cog(Hybrid(bot))
