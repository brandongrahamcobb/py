# Main.py

import asyncio
import logging
import logging.handlers
import os
import yaml

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

    def load_token():
        config_file = 'config.yaml'
        if os.path.exists(config_file):
            with open(config_file, 'r') as file:
                config = yaml.safe_load(file)
                return config.get('token')
        return None

    def save_token(token):
        config = {'token': token}
        with open('config.yaml', 'w') as file:
            yaml.safe_dump(config, file)

    async def validate_token(token):
        intents = discord.Intents.default()
        bot = commands.Bot(command_prefix='!', intents=intents)
        try:
            await bot.login(token)
            await bot.close()
            return True
        except discord.LoginFailure:
            return False

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
        TOKEN = Main.load_token()
        if not TOKEN or not await Main.validate_token(TOKEN):
            while True:
                TOKEN = input("Enter your Discord bot token: ").strip()
                if await Main.validate_token(TOKEN):
                    Main.save_token(TOKEN)
                break
            else:
                print("Invalid token. Please try again.")
        await bot.start(TOKEN)

asyncio.run(main())


