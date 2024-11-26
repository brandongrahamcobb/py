""" Vyrtuous, a discord.py bot, is an open-source package containing three cogs and a helper file.
    Copyright (C) 2024  github.com/brandongrahamcobb

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public license as published by
    the Free Software Foundation, either version 3 of the license, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public license for more details.

    You should have received a copy of the GNU General Public license
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from discord.ext import commands
from typing import List, Optional, Dict, Any
from os import getenv, makedirs
from os.path import abspath, dirname, expanduser, exists, isfile, join
from bot.utils import helpers

import asyncio
import bot.utils.helpers as helpers
import code
import discord
import os
import threading
import sys
import yaml

global logger

home = expanduser('~')

def prompt_for_values(prompt: str, default_value: str) -> str:
    value = input(f'{prompt} [{default_value}]: ')
    return value if value else default_value

class Vyrtuous(commands.Bot):

    _config = None                         # This is to establish a class variable.

    def __init__(
        self,
        *args,
        initial_extensions: List[str],
        testing_guild_id: Optional[int] = None,
        **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.testing_guild_id = testing_guild_id
        self.initial_extensions = initial_extensions

    @classmethod
    def _get_config(cls) -> Dict[str, Any]:
        if cls._config is None:
            if isfile(helpers.path_config_yaml):
                with open(helpers.path_config_yaml, 'r') as file:
                    data = yaml.safe_load(file)
                data['cogs'] = prompt_for_values('Enter the cogs.', data.get('cogs', [
                    'bot.cogs.hybrid',
                    'bot.cogs.indica',
                    'bot.cogs.sativa'
                ]))
                data['command_prefix'] = prompt_for_values('Enter the command prefix', data.get('command_prefix', '!'))
                data['database_url'] = prompt_for_values('Enter the database URL', data.get('database_url', 'brandongcobb.com'))
                data['intents'] = prompt_for_values('Enter the intents', data.get('intents', 'discord.Intents.all()'))
                data['logging_level'] = prompt_for_values('Enter the logging level', data.get('logging_level', 'INFO'))
                data['owner_id'] = prompt_for_values('Enter the owner ID', data.get('owner_id', '154749533429956608'))
                data['testing_guild_id'] = prompt_for_values('Enter the testing guild ID', data.get('testing_guild_id', '1286307083976835103'))
                data['token'] = prompt_for_values('Enter the bot token', data.get('token', ''))
                data['user_agent'] = prompt_for_values('Enter the User-Agent header', data.get('user_agent', 'Vyrtuous'))
                data['version'] = prompt_for_values('Enter the bot version', data.get('version', '1.0.0'))
                data['api_keys'] = data.get('api_keys', {})
                for i in range(1, 21):
                    key = f'api_key_{i}'
                    current_key = data['api_keys'].get(key, '')
                    data['api_keys'][key] = prompt_for_values(f'Enter API key {i}', current_key)
            else:
                makedirs(dirname(helpers.path_config_yaml), exist_ok=True)
                data = {
                    'api_keys': {},
                    'cogs': prompt_for_values('Enter the cogs.', [
                        'bot.cogs.hybrid',
                        'bot.cogs.indica',
                        'bot.cogs.sativa'
                    ]),
                    'command_prefix': prompt_for_values('Enter the command prefix', '!'),
                    'database_url': prompt_for_values('Enter the database URL', 'brandongcobb.com'),
                    'intents': prompt_for_values('Enter the intents', 'discord.Intents.all()'),
                    'logging_level': prompt_for_values('Enter the logging level', 'INFO'),
                    'owner_id': prompt_for_values('Enter the your user ID', '154749533429956608'),
                    'testing_guild_id': prompt_for_values('Enter the testing guild ID', '1286307083976835103'),
                    'token': prompt_for_values('Enter the bot token', ''),
                    'version': prompt_for_values('Enter the bot version', '1.0.0'),
                    'user_agent': prompt_for_values('Enter the User-Agent header', 'Vyrtuous'),
                }
                for i in range(1, 21):
                    data['api_keys'][f'api_key_{i}'] = prompt_for_values(f'Enter API key {i}', '')
            with open(helpers.path_config_yaml, 'w') as file:
                yaml.dump(data, file)
            cls._config = data
            return cls._config

    async def load_config(self) -> Dict[str, Any]:
        return self._get_config()

    async def setup_hook(self) -> None:
        for cog in self.initial_extensions:
            await self.load_extension(cog)
        if self.testing_guild_id:
            guild = discord.Object(self.testing_guild_id)
            self.tree.copy_global_to(guild=guild)
            await self.tree.sync(guild=guild)

async def main():

    intents = discord.Intents.all()

    async with Vyrtuous(
        command_prefix='!',
        intents=intents,
        initial_extensions=[
            'bot.cogs.hybrid',
            'bot.cogs.indica',
            'bot.cogs.sativa'
        ],
    ) as bot:

        config = await bot.load_config()
        bot.config = config

        bot.helpers = helpers
        helpers.increment_version(config)
        helpers.setup_logging(config)
        await bot.start(config['token'])


def run():
    asyncio.run(main())

if __name__ == '__main__':
    run()

