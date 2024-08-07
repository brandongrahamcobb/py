""" main.py
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
"""

from discord.ext import commands
from typing import List, Optional
from utils.helpers import load_config
from utils.helpers import setup_logging

import asyncio
import discord
import os
#import json
#import requests
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

class Lucy(commands.Bot):
    def __init__(
        self,
        *args,
        initial_extensions: List[str],
        testing_guild_id: Optional[int] = None,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.initial_extensions = initial_extensions
        self.testing_guild_id = testing_guild_id
    async def setup_hook(self) -> None:
        for extension in self.initial_extensions:
            await self.load_extension(extension)
        if self.testing_guild_id:
            guild = discord.Object(self.testing_guild_id)
            self.tree.copy_global_to(guild = guild)
            await self.tree.sync(guild = guild)

async def main():
    config = load_config()
    token = config['token']
    setup_logging()
    exts = ['cogs.admin_cog', 'cogs.game_cog', 'cogs.my_cog', 'cogs.user_cog']
    intents = discord.Intents.all()
    intents.message_content = True
    async with Lucy(
        commands.when_mentioned_or('!'),
        initial_extensions = exts,
        intents = intents,
    ) as bot:
        await bot.start(token)

asyncio.run(main())
