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
from bot.utils.helpers import load_config
from bot.utils.helpers import setup_logging

import asyncio
import discord

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
    async with Lucy(
        commands.when_mentioned_or(config['command_prefix']),
        initial_extensions = config.get('cogs', []),
        intents = eval(config['intents']),
        testing_guild_id = config['testing_guild_id']
    ) as bot:
        await bot.start(token)

        @bot.event
        async def on_message(message):
            if message.author.bot:
                return
            # CobbBrandonGraham_helpers_110824.py
            add_user_if_new(str(message.author.id), message.author.name)

if __name__ == "__main__":
    asyncio.run(main())
