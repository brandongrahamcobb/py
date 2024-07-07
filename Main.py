# Main.py

import asyncio
import logging
import logging.handlers
import os

from typing import List, Optional

import discord
import dotenv
from discord.ext import commands

class Client(commands.Bot):
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
    async with Client(
        commands.when_mentioned_or('!'),
        initial_extensions = exts,
        intents = intents,
    ) as bot:
        dotenv.load_dotenv()
        await bot.start(os.getenv("DISCORD_TOKEN"))

asyncio.run(main())


