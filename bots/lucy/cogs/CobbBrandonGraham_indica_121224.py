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
from collections import defaultdict
from discord.ext import commands, tasks
from gtts import gTTS
from os.path import abspath, dirname, expanduser, join

import asyncio
import bot.utils.chem_helpers as chem_helpers
import bot.utils.helpers as helpers
import bot.utils.openai_helpers as openai_helpers
import datetime
import discord
import json
import opuslib
import os
import subprocess
import wave
import traceback

dir_base = dirname(abspath(__file__))
path_chem_helpers = join(dir_base, '..', 'utils', 'chem_helpers.py')
path_openai_helpers = join(dir_base, '..', 'utils', 'openai_helpers.py')
path_home = expanduser('~')
path_hybrid = join(dir_base, '..', 'cogs', 'hybrid.py')
path_indica = join(dir_base, '..', 'cogs', 'indica.py')
path_main = join(dir_base, '..', 'main.py')
path_sativa = join(dir_base, '..', 'cogs', 'sativa.py')

class Indica(commands.Cog):

    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.hybrid = self.bot.get_cog('Hybrid')
        self.sativa = self.bot.get_cog('Sativa')
        self.chem_helpers_py = helpers.load_contents(path_chem_helpers)
        self.openai_helpers_py = helpers.load_contents(path_openai_helpers)
        self.hybrid_py = helpers.load_contents(path_hybrid)
        self.indica_py = helpers.load_contents(path_indica)
        self.main_py = helpers.load_contents(path_main)
        self.sativa_py = helpers.load_contents(path_sativa)
        self.sys_input = f"""
            Your main.py file is {self.main_py}.
            Your cogs are in cogs/ {self.hybrid_py}, {self.indica_py}, {self.sativa_py}.
            Your chem_helpers are in utils/ {self.chem_helpers_py}.
            Your openai_helpers are in utils/ {self.openai_helpers_py}.
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
        if message.author == self.bot.user or message.content.startswith('!'):
            return
        try:
            if message.content:
#                async for moderation in openai_helpers.create_pseudomoderation(input_text=message.content, conversation_id='TextModerator', sys_input=await self.add_application_info("")):
 #                   await openai_helpers.handle_api_moderation(message, moderation)
                return
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
                        async for moderation in openai_helpers.create_moderation(input_text):
                            await openai_helpers.handle_api_moderation(message, moderation)
        except Exception as e:
            print(e)

    @commands.Cog.listener()
    async def on_ready(self):
        bot_user = self.bot.user
        bot_name = bot_user.name
        bot_id = bot_user.id
        guild_count = len(self.bot.guilds)
        info = (
            f'\n=============================\n'
            f'bot Name: {bot_name}\n'
            f'bot ID: {bot_id}\n'
            f'Connected Guilds: {guild_count}\n'
            f'============================='
        )
        guild_info = '\n'.join(
            [f'- {guild.name} (ID: {guild.id})' for guild in self.bot.guilds]
        )
        stats_message = f'{info}\n\nGuilds:\n{guild_info}'
        print(stats_message)
        user = await self.bot.fetch_user(self.config['owner_id'])
        await user.send(stats_message)

async def setup(bot: commands.Bot):
    await bot.add_cog(Indica(bot))
