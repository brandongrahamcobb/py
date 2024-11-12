''' script_cog.py The purpose of this program is to provide religious scripture referencing to a Discord bot.
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
'''

from bs4 import BeautifulSoup
from discord import app_commands, Embed
from discord.ext import commands
from typing import List

import bot.utils.helpers as helpers

import discord
import bs4
import io
import os
import requests
import shlex
import traceback

class ScriptCog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config

    @commands.command(name='script', description='Usage !script <NIV/ESV> <Book>.<Chapter>.<Verse>')
    async def script(self, ctx: commands.Context, version: str, *, reference: str):
         try:
             BIBLE_IDS = {
                 'esv': 'de4e12af7f28f599-02',
                 'nkjv': 'de4e12af7f28f599-01',
                 'niv': '06125adad2d5898a-01',
             }
             version = version.lower()
             if version in BIBLE_IDS:
                 bible_id = BIBLE_IDS[version]
                 api = f'https://api.scripture.api.bible/v1/bibles/{bible_id}/search?query={reference}'
                 response = requests.get(api, headers=helpers.get_scripture_headers())
                 if response.ok:
                     json = response.json()
                     passages = json.get('data', {}).get('passages', [])
                     soup = BeautifulSoup(passages[0].get('content'), 'html.parser')
                     soup.get_text()
                     cleaned_content = soup.get_text()
                     message = f'**{reference}** ({version.upper()})\n{cleaned_content}'
                     await ctx.send(message)
             if version == '':
                 response = requests.get(f'https://api.alquran.cloud/v1/ayah/{reference}/en.asad', headers=get_scripture_headers())
                 if response.ok:
                     json = response.json()
                     message = f'**{reference}** ({version.upper()})\n{json['data']['text']}'
                     await ctx.send(message)
         except Exception as e:
            await Command.warn(self, ctx, e)

async def setup(bot: commands.Bot):
    await bot.add_cog(ScriptCog(bot))
