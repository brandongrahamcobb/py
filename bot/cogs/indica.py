''' ListenerCog.py The purpose of the program is to contain a declaration of a few self variables and all the listeners for a discord bot.
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
from discord.ext import commands, tasks
from os.path import abspath, dirname, expanduser, join

import asyncio
import datetime
import discord
import os
import subprocess
import traceback

CHECKMARK_EMOJI = '✅'
CROSS_EMOJI = '❌'
ROLE_ID = 1308689505158565918
dir_base = dirname(abspath(__file__))
path_openai_helpers = join(dir_base, '..', 'utils', 'openai_helpers.py')
path_home = expanduser('~')
path_main = join(dir_base, '..', 'main.py')
path_sativa = join(dir_base, '..', 'cogs', 'sativa.py')

class Indica(commands.Cog):

    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.hybrid = self.bot.get_cog('Hybrid')
        self.sativa = self.bot.get_cog('Sativa')
        self.add_watermark = load_contents(helpers.PATH_ADD_WATERMARK)
        self.adjust_hue_and_saturation = load_contents(helpers.PATH_ADJUST_HUE_AND_SATURATION)
        self.arpp = load_contents(helpers.PATH_ARPP)
        self.combine = load_contents(helpers.PATH_COMBINE)
        self.config = load_contents(helpers.PATH_CONFIG)
        self.create_batch_completion = load_contents(helpers.PATH_CREATE_BATCH_COMPLETION)
        self.create_completion_deprecated = load_contents(helpers.PATH_CREATE_COMPLETION_DEPRECATED)
        self.create_completion = load_contents(helpers.PATH_CREATE_COMPLETION)
        self.create_https_completion = load_contents(helpers.PATH_CREATE_HTTPS_COMPLETION)
        self.create_moderation = load_contents(helpers.PATH_CREATE_MODERATION)
        self.discord = load_contents(helpers.PATH_DISCORD)
        self.draw_fingerprint = load_contents(helpers.PATH_DRAW_FINGERPRINT)
        self.draw_watermarked_molecule = load_contents(helpers.PATH_DRAW_WATERMARKED_MOLECULE)
        self.format_error_check = load_contents(helpers.PATH_FORMAT_ERROR_CHECK)
        self.get_molecule_name = load_contents(helpers.PATH_GET_MOLECULE_NAME)
        self.get_mol = load_contents(helpers.PATH_GET_MOL)
        self.get_proximity = load_contents(helpers.PATH_GET_PROXIMITY)
        self.get_scripture = load_contents(helpers.PATH_GET_SCRIPTURE)
        self.google = load_contents(helpers.PATH_GOOGLE)
        self.gsrs = load_contents(helpers.PATH_GSRS)
        self.helpers = load_contents(helpers.PATH_HELPERS)
        self.hybrid = join(dir_base, '..', 'cogs', 'hybrid.py')
        self.indica = join(dir_base, '..', 'cogs', 'indica.py')
        self.load_contents = load_contents(helpers.PATH_LOAD_CONTENTS)
        self.load_yaml = load_contents(helpers.PATH_LOAD_YAML)
        self.setup_logging = load_contents(helpers.PATH_SETUP_LOGGING)
        self.stable_cascade = load_contents(helpers.PATH_STABLE_CASCADE)
        self.unique_pairs = load_contents(helpers.PATH_UNIQUE_PAIRS)
        self.sum_of_paths = """
            self.adjust_hue_and_saturation + self.arpp + self.bot + self.combine +
            self.create_batch_completion + self.create_completion +
            self.create_https_completion + self.create_moderation +
            self.draw_fingerprint + self.draw_watermarked_molecule +
            self.format_error_check + self.get_molecule_name + self.get_mol +
            self.get_proximity + self.get_scripture + self.google + self.gsrs +
            self.load_contents + self.load_yaml + self.setup_logging +
            self.stable_cascade + self.unique_pairs + self.vyrtuous
        """
        self.sys_input = f"""
            Your utilities are {self.sum_of_paths}.
        """
    
    async def add_application_info(self, text_input):
        text = await self.bot.application_info()
        text_full = text.description + text_input
        return text_full

    @commands.Cog.listener()
    async def on_message_edit(self, before, after):
        if before.content != after.content:
            ctx = await self.bot.get_context(after)
            if ctx.command:
                await self.bot.invoke(ctx)

    @commands.Cog.listener()
    async def on_message(self, message):
        try:
            if message.attachments:
                for attachment in message.attachments:
                    if attachment.content_type.startswith('image/'):
                        input_text = [
                            {
                                "type": "text",
                                "text": message.content
                            },
                            {
                                "type": "image_url",
                                "image_url": {
                                    'url': attachment.url
                                },
                            },
                        ]
                        async for response in openai_helpers.create_completion(completions=1, conversation_id='ChemistryWizard', input_text=input_text, max_tokens=16384, model='gpt-4o-mini', stop='', store=False, stream=False, sys_input=await self.add_application_info(self.sys_input), temperature=1.0, top_p=1.0):
                            await message.channel.send(response)
                        #async for moderation in openai_helpers.create_moderation(input_text):
                         #   await openai_helpers.handle_api_moderation(message, moderation)
        except Exception as e:
            print(e)

   @commands.Cog.listener()
    async def on_ready(self):
        bot_user = self.bot.user
        bot_name = bot_user.name
        bot_id = bot_user.id
        guild_count = len(self.bot.guilds)
        embed = discord.Embed(
            title="Bot Status",
            description="The bot is now online and ready!",
            color=discord.Color.green()
        )
        embed.add_field(name="Bot Name", value=f"{bot_name}", inline=False)
        embed.add_field(name="Bot ID", value=f"{bot_id}", inline=False)
        embed.add_field(name="Connected Guilds", value=f"{guild_count}", inline=False)
        guild_info = '\n'.join([f"- {guild.name} (ID: {guild.id})" for guild in self.bot.guilds])
        embed.add_field(name="Guilds", value=guild_info if guild_info else "No guilds connected", inline=False)
        channel_id = 1315735859848544378
        channel = self.bot.get_channel(channel_id)
        if channel:
            await channel.send(embed=embed)
        print("=============================")
        print(f"Bot Name: {bot_name}")
        print(f"Bot ID: {bot_id}")
        print(f"Connected Guilds: {guild_count}")
        print("Guilds:")
        print(guild_info)
        print("=============================")
        owner_id = self.bot.config.get('owner_id')  # Ensure you have 'owner_id' in config
        if owner_id:
            try:
                owner = await self.bot.fetch_user(owner_id)
                await owner.send(embed=embed)
            except Exception as e:
                print(f"Could not DM owner: {e}")
async def setup(bot: commands.Bot):
    await bot.add_cog(Indica(bot))
