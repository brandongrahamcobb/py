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
import bot.utils.helpers as helpers
import datetime
import discord
import opuslib
import os
import subprocess
import wave

dir_base = dirname(abspath(__file__))
path_helpers = join(dir_base, '..', 'utils', 'helpers.py')
path_home = expanduser('~')
path_hybrid = join(dir_base, '..', 'cogs', 'hybrid.py')
path_indica = join(dir_base, '..', 'cogs', 'indica.py')
path_main = join(dir_base, '..', 'main.py')
path_sativa = join(dir_base, '..', 'cogs', 'sativa.py')
path_users_yaml = join(path_home, '.config', 'vyrtuous', 'users.yaml')

import json

class Indica(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.hybrid = self.bot.get_cog('Hybrid')
        self.sativa = self.bot.get_cog('Sativa')
        self.helpers_py = helpers.load_contents(path_helpers)
        self.hybrid_py = helpers.load_contents(path_hybrid)
        self.indica_py = helpers.load_contents(path_indica)
        self.main_py = helpers.load_contents(path_main)
        self.sativa_py = helpers.load_contents(path_sativa)
        self.hybrid = self.bot.get_cog('Hybrid')
        self.sativa = self.bot.get_cog('Sativa')
        self.sys_input = f"""
            Your main.py file is {self.main_py}.
            Your cogs are in cogs/ {self.hybrid_py}, {self.indica_py}, {self.sativa_py}.
            Your helpers are in utils/ {self.helpers_py}.
        """
        self.channel_id = 1305608084017905734
#        self.post_message.start()
        self.users = {}
    
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
        # Ignore messages from self.
        if message.author == self.bot.user or message.content.startswith('!'):
            return
        if message.author == self.bot.user or '!' in message.content[0]:
            return
        if self.bot.user.mentioned_in(message) or isinstance(message.channel, discord.DMChannel):
            if int(message.guild.id) == int(self.config['testing_guild_id']):
                conversation_id = message.channel.id + message.author.id
                async for response in helpers.deprecated_create_completion(f'{message.content}', await self.add_application_info(self.sys_input), conversation_id):
                    responses = helpers.chunk_string(text=response)
                    for response in responses:
                        await message.channel.send(f"@{message.author.name}, {response}")
        if int(message.guild.id) == int(self.config['testing_guild_id']):
            async for moderation in helpers.create_moderation(message.content):
                if isinstance(moderation, dict) and 'error' in moderation:
                    print(f"Moderation error: {moderation['error']}")
                    return  # Handle errors (e.g., log or notify users)
            # Check the structure of the moderation response
                if 'results' in moderation and len(moderation['results']) > 0:
                    flagged = moderation['results'][0].get('flagged', False)
                    print("Flagged Status:", flagged)
                    categories = moderation['results'][0].get('categories', {})
                    print("Categories:", categories)
                    if flagged:
                        await message.delete()  # Delete the flagged message
                        await message.channel.send("Your message was flagged and deleted due to inappropriate content.")
                else:
                    print("No moderation results found.")
        if int(message.guild.id) == int(self.config['testing_guild_id']):
            moderation_input = f'{message.content} in channel: {message.channel.name}'
            structured_sys_input = self.get_sys_input(message.id)
            try:
                async for moderation_response in helpers.deprecated_create_moderation(
                    input_text=moderation_input,
                    sys_input=structured_sys_input,
                    conversation_id='1'
                ):
                    if not moderation_response or not moderation_response.strip():
                        print("Received an empty or null response from moderation.")
                        continue  # Skip this iteration and continue with the next
                    try:
                        moderation_dict = json.loads(moderation_response)
                    except json.JSONDecodeError as e:
                        print(f"JSON decoding failed: {e}")
                        print(f"Response content was: '{moderation_response}'")
                        continue  # Skip this iteration and continue with the next
                    flagged = moderation_dict.get('results', [{}])[0].get('flagged', False)
                    await self.handle_moderation_result(flagged, message)
            except Exception as e:
                print(f"Error during moderation: {e}")

    def get_sys_input(self, message_id):
        # Provide clear instructions for both input and expected output
        return f"""
You are a moderation assistant. Please analyze the following message and respond in this structured JSON format:
{{
    "id": "{message_id}",
    "model": "gpt-4o-mini",
    "results": [
        {{
            "flagged": false,
            "categories": {{
                "sexual": false,
                "sexual/minors": false,
                "harassment": false,
                "harassment/threatening": false,
                "violence": false,
                "violence/graphic": false,
                "self-harm": false,
                "self-harm/intent": false,
                "self-harm/instructions": false,
                "self-harm/ideation": false,
                "hate": false,
                "hate/threatening": false,
                "extremism": false,
                "extremism/violence": false
            }},
            "category_scores": {{
                "sexual": 0,
                "sexual/minors": 0,
                "harassment": 0,
                "harassment/threatening": 0,
                "violence": 0,
                "violence/graphic": 0,
                "self-harm": 0,
                "self-harm/intent": 0,
                "self-harm/instructions": 0,
                "self-harm/ideation": 0,
                "hate": 0,
                "hate/threatening": 0,
                "extremism": 0,
                "extremism/violence": 0
            }}
        }}
    ]
}}
"""

    async def handle_moderation_result(self, flagged, message):
        if flagged:
            await message.delete()
            count = helpers.increment_infraction(str(message.author.id), self.users)
            await helpers.save_yaml(path_users_yaml, self.users)
            embed = helpers.create_embed(
                "Message Deleted",
                f"Your message was deleted due to inappropriate content. Infractions: {count}"
            )
            await message.author.send(embed=embed)
            # This may also call for a response
#            async for response in helpers.deprecated_create_completion(f'{message.content}', self.sys_input, conversation_id=None):
#                responses = helpers.chunk_string(text=response)
#                await message.reply(f"@{message.author.name}, {responses[0]}")
#                for response in responses[1:]:
#                    await message.reply(response)

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
