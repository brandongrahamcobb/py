# Main.py

import asyncio
import logging
import logging.handlers
import os
import json
import requests

from typing import List, Optional

import discord
from discord.ext import commands

class Main(commands.Bot):
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

    def verify_discord_token(token):
        headers = {
           'Authorization': f'Bot {token}'
        }
        url = 'https://discord.com/api/v9/users/@me'

        try:
            response = requests.get(url, headers=headers)
            response.raise_for_status()  # Raise an exception for 4xx or 5xx errors
            data = response.json()
            return data
        except requests.exceptions.HTTPError as err:
            print(f"HTTP error occurred: {err}")
            return None
        except Exception as err:
            print(f"Error occurred: {err}")
            return None

    def read_token_from_file(filename):
        if os.path.exists(filename):
            with open(filename, 'r') as f:
                data = json.load(f)
                return data.get('token', None)
        return None

    def write_token_to_file(filename, token):
        data = {'token': token}
        with open(filename, 'w') as f:
            json.dump(data, f)

async def main():
    logger = logging.getLogger('discord')
    logger.setLevel(logging.INFO)

    handler = logging.handlers.RotatingFileHandler(
       filename='../log/discord.log',
       encoding='utf-8',
       maxBytes=32 * 1024 * 1024,
       backupCount=5,
    )
    dt_fmt = '%Y-%m-%d %H:%M:%S'
    formatter = logging.Formatter('[{asctime}] [{levelname:<8} {name}: {message}', dt_fmt, style='{')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    exts = ['Chemistry', 'Listener', 'Super']
    intents = discord.Intents.all()
    intents.message_content = True
    async with Main(
        commands.when_mentioned_or('!'),
        initial_extensions = exts,
        intents = intents,
    ) as bot:
        token_filename = 'config.json'
        token = Main.read_token_from_file(token_filename)
        while token is None or len(token) == 0 or Main.verify_discord_token(token) is None:
            token = input("Enter your Discord bot token: ").strip()
            Main.write_token_to_file(token_filename, token)
            print(f"Token saved to {token_filename}")
        await bot.start(token)

asyncio.run(main())


