"""
    ai_cog.py The purpose of this program is to provide functionality to a nicotine replacement therapy bot on Discord
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
"""

from collections import defaultdict
from discord.ext import commands
from openai import AsyncOpenAI

import asyncio
import discord
import traceback

class GPT:
    def __init__(self, api_key):
        self.ai_client = AsyncOpenAI(api_key=api_key)
        self.conversations = defaultdict(list)

    async def ai(self, input_text, sys_input, conversation_id):
        try:
            messages = self.conversations[conversation_id]
            messages.append({'role': 'system', 'content': sys_input})
            messages.append({'role': 'user', 'content': input_text})
            stream = await self.ai_client.chat.completions.create(
                model='gpt-4o-mini',
                messages=messages,
                stream=True
            )
            full_response = ''
            async for chunk in stream:
                content = chunk.choices[0].delta.content
                if content is not None:
                    full_response += content
            self.conversations[conversation_id].append({'role': 'assistant', 'content': full_response})
            yield full_response
        except Exception as e:
            yield traceback.format_exc()

class Researcher(commands.Cog):

    def __init__(self, bot):
       self.bot = bot
       self.config = bot.config
       self.ai = GPT(self.config['api_keys']['api_key_2'])

    @commands.Cog.listener()
    async def on_message(self, message):
        conversation_id = message.channel.id
        if message.author == self.bot.user:
            return
        if self.bot.user.mentioned_in(message) or isinstance(message.channel, discord.DMChannel):
            async with message.channel.typing():
                async for response in self.ai.ai(f'{message.content}', 'I want you to be biased. Your bias is a presciption of veganism. Respond with an empty string if little value is provided by responding. Respond with less than 1500 characters.', conversation_id):
                    await message.channel.send(response)

async def setup(bot):
    await bot.add_cog(Researcher(bot))
