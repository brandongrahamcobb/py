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
from typing import List

import asyncio
import discord

class CustomBot(commands.Bot):
    def __init__(
        self,
        *args,
        initial_extension: List[str],
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.initial_extension = initial_extension

    async def setup_hook(self) -> None:
        for cog in self.initial_extension:
            await self.load_extension(cog)

async def main():
    token = YOUR_TOKEN_HERE
    intents = discord.Intents.all()
    intents.message_content = True
    initial_extension = 'cogs.my_cog',
    async with CustomBot(
        commands.when_mentioned_or('\\'),
        initial_extension=initial_extension,
        intents=intents,
    ) as bot:
        await bot.start(token)

def run():
    asyncio.run(main())

if __name__ == '__main__':
    run()
