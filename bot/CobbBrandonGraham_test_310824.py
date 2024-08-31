""" main.py
    Copyright (C) 2024 github.com/<put-your-github-here>

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
from os import makedirs
from os.path import dirname, expanduser, isfile, join
from typing import Any, Dict, List

import asyncio
import discord
import os
import yaml

def prompt_for_values(prompt: str, default_value: str) -> str:
    value = input(f'{prompt} [{default_value}]: ')
    return value if value else default_value

class CustomBot(commands.Bot):

    _config = None

    def __init__(
        self,
        *args,
        initial_extension: List[str],
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.initial_extension = initial_extension

    @classmethod
    def _get_config(cls) -> Dict[str, Any]:
        if cls._config is None:
            configuration = join(expanduser('~'), '.config', 'custom_bot', 'config')
            if isfile(configuration):
                with open(configuration, 'r') as file:
                    data = yaml.safe_load(file)
                data['command_prefix'] = prompt_for_values('Please press enter to keep this command prefix or enter a new one.', data.get('command_prefix', '!'))
                data['owner_id'] = prompt_for_values('Please press enter to keep this owner id or enter a new one.', data.get('owner_id', ''))
                data['testing_guild_id'] = prompt_for_values('Please press enter to keep this testing guild id or enter a new one.', data.get('testing_guild_id', '1110691275948183634'))
                data['token'] = prompt_for_values('Please press enter to keep this token or enter a new one.', data.get('token', ''))
                data['api_keys'] = data.get('api_keys', {})
                for i in range(1, 21):
                    key = f'api_key_{i}'
                    current_key = data['api_keys'].get(key, '')
                    data['api_keys'][key] = prompt_for_values(f'Enter API key {i}', current_key)
            else:
                makedirs(dirname(configuration), exist_ok=True)
                data = {
                    'api_keys': {},
                    'command_prefix': prompt_for_values('Enter the command prefix, or press enter to skip', '!'),
                    'owner_id': prompt_for_values('Enter your user ID, or press enter to skip.', ''),
                    'testing_guild_id': prompt_for_values('Enter the testing guild ID, or press enter to skip.', '1110691275948183634'),
                    'token': prompt_for_values('Enter the bot token.', ''),
                }
                for i in range(1, 21):
                    data['api_keys'][f'api_key_{i}'] = prompt_for_values(f'Enter API key {i}, or press enter to skip.', '')
            with open(configuration, 'w') as file:
                yaml.dump(data, file)
            cls._config = data
            return cls._config

    async def load_config(self) -> Dict[str, Any]:
        return self._get_config()

    async def setup_hook(self) -> None:
        for cog in self.initial_extension:
            await self.load_extension(cog)

async def main():

    async with CustomBot(
        command_prefix=commands.when_mentioned_or('!'),
        intents=discord.Intents.all(),
        initial_extension=[
            'bot.cogs.my_cog'
        ],
    ) as bot:

        config = await bot.load_config()
        bot.config = config

        await bot.start(config['token'])

def run():
    asyncio.run(main())

if __name__ == '__main__':
    run()
