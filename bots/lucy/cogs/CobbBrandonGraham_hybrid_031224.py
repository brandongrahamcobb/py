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

import asyncio
import discord
import io
import os
import shlex
import traceback

class Hybrid(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.indica = self.bot.get_cog('Indica')
        self.sativa = self.bot.get_cog('Sativa')

        self.sys_input = """
         """

    @commands.command(description='Change your role color using RGB values. Usage: between `!colorize 0 0 0` and `!colorize 255 255 255`')
    async def colorize(self, ctx: commands.Context, r: int = commands.parameter(default="149", description="Anything between 0 and 255."), g: int = commands.parameter(default="165", description="Anything betwen 0 and 255."), b: int = commands.parameter(default="165", description="Anything between 0 and 255.")):
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

    @commands.command(name='colors')
    async def colors(self, ctx: commands.Context, *args):
        try:
            attachment = ctx.message.attachments[0]
            img_data = requests.get(attachment.url).content
            image = Image.open(io.BytesIO(img_data))
            image = image.resize((100, 100))
            image = image.convert('RGB')
            pixels = list(image.getdata())
            pixels = [pixel for pixel in pixels if not (pixel[0] > 150 and pixel[1] > 150 and pixel[2] > 150)]
            pixels = [pixel for pixel in pixels if not (pixel[0] < 10 and pixel[1] < 10 and pixel[2] < 10)]
            color_counts = Counter(pixels)
            predominant_colors = color_counts.most_common(int(args[0]))
            message = "Predominant colors (excluding whites above (150, 150, 150)):\n"
            for color, count in predominant_colors:
                message += f"Color: {color} Count: {count}\n"
            await ctx.send(message)
        except Exception as e:
            await ctx.send(e)

    @commands.hybrid_command(name='get', description='Limited usage. Usage: !get <prompt-for-image>.')
    async def get(self, ctx: commands.Context, *, argument: str = commands.parameter(default=None, description="Image prompt.")):
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

    @commands.hybrid_command(name='join', description='Join a voice channel. Usage: !join [channel_id]')
    async def join(self, ctx: commands.Context, channel_id: int):
        try:
            voice_client = await helpers.join_voice_channel(channel_id)
            await ctx.send(f"Joined the voice channel: {voice_client.channel.name}")
        except Exception as e:
            await ctx.send(str(e))

    @commands.hybrid_command(name='leave', description='Leave the voice channel. Usage: !leave')
    async def leave(self, ctx: commands.Context):
        voice_client = discord.utils.get(self.bot.voice_clients)
        if voice_client and voice_client.is_connected():
            await voice_client.disconnect()
            await ctx.send("Disconnected from the voice channel.")
        else:
            await ctx.send("I'm not currently in a voice channel.")

    @commands.command(name='load', hidden=True)
    async def load(self, ctx: commands.Context, *, module: str):
        try:
            await ctx.bot.load_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')

    @commands.hybrid_command(name='draw', description='Usage: !draw glow <molecule> or !draw gsrs <molecule> or !draw shadow <molecule>. Use quotations for multistring molecules.')
    async def molecule(self, ctx: commands.Context, option: str = commands.parameter(default="glow", description="Compare `compare or Draw style `glow` `gsrs` `shadow`."), *, molecules: str = commands.parameter(default=None, description="Any molecule"), quantity: int = commands.parameter(default=1, description="Quantity of glows")):
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
                for _ in range(quantity.default):
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

    @commands.command(name='pomodoro', help='Start a Pomodoro timer.')
    async def pomodoro(self, ctx, work_time: int = 25, short_break: int = 5, long_break: int = 15, cycles: int = 4):
        await ctx.send(f"Starting a Pomodoro timer with {work_time} minutes work, {short_break} minutes short break, "
                       f"and {long_break} minutes long break for {cycles} cycles.")
        for _ in range(cycles):
            await ctx.send(f"**Work Session**: {work_time} minutes!")
            await asyncio.sleep(work_time * 60)
            await ctx.send(f"**Short Break**: {short_break} minutes!")
            await asyncio.sleep(short_break * 60)
        await ctx.send(f"**Long Break**: {long_break} minutes!")
        await asyncio.sleep(long_break * 60)
        await ctx.send("Pomodoro session completed. Great job!")

    @commands.hybrid_command(name='talk', description='Start talking. Usage: !talk')
    async def talk(self, ctx: commands.Context):
        voice_client = discord.utils.get(self.bot.voice_clients)
        if voice_client:
            await ctx.send("I am ready to listen for 30 seconds!")
            user_input = await helpers.voice_to_text(ctx)

            if user_input:
                async for response in helpers.deprecated_create_completion(user_input, self.sys_input, 'VC'):
                    if response:
                        await ctx.send(response)  # Send the response message
                        await helpers.text_to_voice(response, voice_client)  # Convert to voice and play it
                else:
                    await ctx.send("I received no response.")
        else:
            await ctx.send("I'm not connected to a voice channel.")

    @commands.command(name='recipe', description='Usage: !recipe')
    async def recipe(self, ctx: commands.Context, *args):
        try:
            if ctx.interaction:
                await ctx.interaction.response.defer(ephemeral=True)
            embed = helpers.get_recipe(args[0])
            await ctx.send(embed=embed)
        except Exception as e:
            await ctx.send(f'An error occurred: {e}')

    @commands.command(name='recipe2', description='Usage: !recipe2')
    async def recipe_two(self, ctx: commands.Context, *args):
        try:
            if ctx.interaction:
                await ctx.interaction.response.defer(ephemeral=True)
            embed = await helpers.get_other_recipe(args[0])
            await ctx.send(embed=embed)
        except Exception as e:
            await ctx.send(f'An error occurred: {e}')

    @commands.hybrid_command(hidden=True)
    async def reload(self, ctx: commands.Context, *, module: str):
        try:
            if ctx.interaction:
                await ctx.interaction.response.defer(ephemeral=True)
            await ctx.bot.reload_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')

    @commands.hybrid_command(name='script', description='Usage !script <ESV/NIV/NKJV> <Book>.<Chapter>.<Verse>. !script Quran <Al-Surah>')
    async def script(self, ctx: commands.Context, version: str = commands.parameter(default="esv", description="ESV/NIV/NKJV/Quran"), *, reference: str = commands.parameter(default=None, description=">Book>.<Chapter>.<Verse> for Bible and <Al-Surah Number> for Quran.")):
        if ctx.interaction:
            await ctx.interaction.response.defer(ephemeral=True)
        await ctx.send(helpers.get_scripture(version, reference))

    @commands.hybrid_command(name='search', description='Usage: !search <query>. Search Google.')
    async def search(self, ctx: commands.Context, *, query: str = commands.parameter(default=None, description="Google search a query.")):
        if ctx.interaction:
            await ctx.interaction.response.defer(ephemeral=True)
        results = helpers.scrape_google_search(query)
        embed = discord.Embed(title=f"Search Results for '{query}'", color=discord.Color.blue())
        for result in results:
            embed.add_field(name=result["title"], value=result["link"], inline=False)
        await ctx.send(embed=embed)

    @commands.hybrid_command(name='warn', description='Usage: !warn. Respond to any message with !warn.')
    async def warn(self, ctx: commands.Context):
        try:
            if ctx.interaction:
                await ctx.interaction.response.defer(ephemeral=True)
            if ctx.message.reference is None:
                await ctx.send("Please reply to a message you want to delete.")
                return
            original_message_id = ctx.message.reference.message_id
            original_message = await ctx.channel.fetch_message(original_message_id)
            try:
                await original_message.delete()  # Attempt to delete the original message
                await ctx.send("Message deleted successfully.")
            except discord.NotFound:
                await ctx.send("The original message was not found.")
            except discord.Forbidden:
                await ctx.send("I do not have permission to delete that message.")
            except discord.HTTPException:
                await ctx.send("An error occurred while trying to delete the message.")
        except Exception as e:
            await ctx.send(f'An error occurred: {e}')

    @commands.command(name='vegan', description='Usage: !vegan')
    async def vegan(self, ctx: commands.Context, *args):
        try:
            if ctx.interaction:
                await ctx.interaction.response.defer(ephemeral=True)
            await ctx.send('https://linktr.ee/veganresource')
        except:
            await ctx.send('Fail')

async def setup(bot: commands.Bot):
    await bot.add_cog(Hybrid(bot))
