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

import asyncio
import logging
import logging.handlers
import os
import json
import requests
import setup

from typing import List, Optional

import discord
from discord.ext import commands

def get_version():
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

def load_config():
    config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'json', 'config.json')
    if os.path.exists(config_path):
        with open(config_path, 'r') as f:
            return json.load(f)
    else:
        raise FileNotFoundError("Configuration file not found.")

# Set up logging
def setup_logging():
    global logger
    log_file = '../log/discord.log'

    # Create a file handler for logging
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)

    # Create a logging format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)

    # Add the file handler to the root logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.addHandler(file_handler)

    # Create the file if it doesn't exist
    if not os.path.exists(log_file):
        open(log_file, 'a').close()

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
    exts = ['game_cog', 'my_cog']
    intents = discord.Intents.all()
    intents.message_content = True
    async with Lucy(
        commands.when_mentioned_or('!'),
        initial_extensions = exts,
        intents = intents,
    ) as bot:
        await bot.start(token)

asyncio.run(main())
