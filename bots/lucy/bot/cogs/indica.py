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
from os.path import abspath, dirname, exists, expanduser, join
from utils.create_https_completion import create_https_completion
from utils.create_moderation import create_moderation
from utils.fine_tuning import TrainingFileBuilder
from utils.load_contents import load_contents
import asyncio
import datetime
import discord
import json
import os
import subprocess
import traceback
import utils.helpers as helpers


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
        self.conversations = defaultdict(list)
        self.temp_conversations = defaultdict(list)
        self.lock = asyncio.Lock()
        self.hybrid = self.bot.get_cog('Hybrid')
        self.sativa = self.bot.get_cog('Sativa')
        self.add_watermark = load_contents(helpers.PATH_ADD_WATERMARK)
        self.adjust_hue_and_saturation = load_contents(helpers.PATH_ADJUST_HUE_AND_SATURATION)
        self.arpp = load_contents(helpers.PATH_ARPP)
        self.benchmark = load_contents(helpers.PATH_BENCHMARK)
        self.clear_screen = load_contents(helpers.PATH_CLEAR_SCREEN)
        self.combine = load_contents(helpers.PATH_COMBINE)
        self.create_batch_completion = load_contents(helpers.PATH_CREATE_BATCH_COMPLETION)
        self.create_completion_deprecated = load_contents(helpers.PATH_CREATE_COMPLETION_DEPRECATED)
        self.create_completion = load_contents(helpers.PATH_CREATE_COMPLETION)
        self.create_https_completion = load_contents(helpers.PATH_CREATE_HTTPS_COMPLETION)
        self.create_moderation = load_contents(helpers.PATH_CREATE_MODERATION)
        self.discord = load_contents(helpers.PATH_DISCORD)
        self.draw_fingerprint = load_contents(helpers.PATH_DRAW_FINGERPRINT)
        self.draw_watermarked_molecule = load_contents(helpers.PATH_DRAW_WATERMARKED_MOLECULE)
        self.fine_tuning = load_contents(helpers.PATH_FINE_TUNING)
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
        self.increment_version = load_contents(helpers.PATH_INCREMENT_VERSION)
        self.load_contents = load_contents(helpers.PATH_LOAD_CONTENTS)
        self.load_yaml = load_contents(helpers.PATH_LOAD_YAML)
        self.prompt_for_values = load_contents(helpers.PATH_PROMPT_FOR_VALUES)
        self.script = load_contents(helpers.PATH_SCRIPT)
        self.setup_logging = load_contents(helpers.PATH_SETUP_LOGGING)
        self.tag = load_contents(helpers.PATH_TAG)
        self.unique_pairs = load_contents(helpers.PATH_UNIQUE_PAIRS)
        self.sum_of_paths = f"""
            {self.adjust_hue_and_saturation} and {self.arpp} and {self.benchmark} and {self.bot} and {self.clear_screen} and {self.combine} and {self.create_batch_completion} and {self.create_completion} and {self.create_https_completion} and {self.create_moderation} and {self.discord} and {self.draw_fingerprint} and {self.draw_watermarked_molecule} and {self.fine_tuning} and {self.format_error_check} and {self.get_molecule_name} and {self.get_mol} and {self.get_proximity} and {self.get_scripture} and {self.google} and {self.gsrs} and {self.helpers} and {self.increment_version} and {self.load_contents} and {self.load_yaml} and {self.setup_logging} and {self.tag} and {self.unique_pairs}
        """
        self.sys_input = f"""
            Your utilities are {self.sum_of_paths}.
        """
    

    @commands.Cog.listener()
    async def on_message_edit(self, before, after):
        if before.content != after.content:
            ctx = await self.bot.get_context(after)
            if ctx.command:
                await self.bot.invoke(ctx)

    @commands.Cog.listener()
    async def on_message(self, message):
        try:
            if self.bot.user == message.author:
                return
            async with self.lock:
                # Input Image/Text
                array = []
                input_text_dict = {
                    "type": "text",
                    "text": message.content.replace('<@1315609784216719370> ', '')
                }
                array.append(input_text_dict)
                if message.attachments:
                    for attachment in message.attachments:
                        if attachment.content_type and attachment.content_type.startswith('image/'):
                            input_image_dict = {
                                'type': 'image_url',
                                'image_url': {
                                    'url': attachment.url
                                }
                            }
                        array.append(input_image_dict)
                    # Image Moderation
                    async for moderation in create_moderation(array):
                        results = moderation.get('results', [])
                        if results and results[0].get('flagged', False):
                            await message.delete()
                            channel = await message.author.create_dm()
                            await channel.send(self.config['openai_moderation_warning'])
                            break
                # Chat Completion
                if self.bot.user in message.raw_mentions and not isinstance(message.type, MessageType.reply):
                    self.conversations.clear()
                if self.bot.user in message.mentions:
                    async for completion in create_https_completion(
                        completions=self.config['openai_chat_n'],
                        conversations=self.conversations,
                        custom_id=message.author.id,
                        input_array=array,
                        max_tokens=self.config['openai_chat_max_tokens'],
                        model=self.config['openai_chat_model'],
                        response_format=self.config['openai_chat_response_format'],
                        stop=self.config['openai_chat_stop'],
                        store=self.config['openai_chat_store'],
                        stream=self.config['openai_chat_stream'],
                        sys_input=self.config['openai_chat_sys_input'], # self.sys_input, 
                        temperature=self.config['openai_chat_temperature'],
                        top_p=self.config['openai_chat_top_p']
                    ):
                        response = completion['choices'][0]['message']['content']
                        self.conversations[message.author.id].append({'role': 'assistant', 'content': response})
                        if len(response) > self.config['discord_character_limit']:
                            await message.reply('My reply was longer than Discord\'s minimum. Oops!')
                        else:
                            await message.reply(response)
                # Chat Moderation
                if self.config['openai_chat_moderation']:
                    async for completion in create_https_completion(
                        completions=helpers.OPENAI_CHAT_MODERATION_N,
                        conversations=self.temp_conversations,
                        custom_id=message.author.id,
                        input_array=input_text_dict,
                        max_tokens=helpers.OPENAI_CHAT_MODERATION_MAX_TOKENS,
                        model=self.config['openai_chat_moderation_model'],
                        response_format=helpers.OPENAI_CHAT_MODERATION_RESPONSE_FORMAT,
                        stop=helpers.OPENAI_CHAT_MODERATION_STOP,
                        store=helpers.OPENAI_CHAT_MODERATION_STORE,
                        stream=helpers.OPENAI_CHAT_MODERATION_STREAM,
                        sys_input=helpers.OPENAI_CHAT_MODERATION_SYS_INPUT,
                        temperature=helpers.OPENAI_CHAT_MODERATION_TEMPERATURE,
                        top_p=helpers.OPENAI_CHAT_MODERATION_TOP_P
                    ):
                        content = json.loads(completion['choices'][0]['message']['content'])
                        flagged = content['results'][0]['flagged']
                        if flagged:
                             channel = await message.author.create_dm()
                             await channel.send(self.config['openai_moderation_warning'])
                             await message.delete()
                # Fine-tuning
        except Exception as e:
            print(e)

async def setup(bot: commands.Bot):
    await bot.add_cog(Indica(bot))
