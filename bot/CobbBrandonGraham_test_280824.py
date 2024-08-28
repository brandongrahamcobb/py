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

import asyncio
import discord

class Lucy(commands.Bot):
    def __init__(
        self,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)

async def main():
    token = YOUR_BOT_TOKEN
    intents = discord.Intents.all()
    intents.message_content = True
    async with Lucy(
        commands.when_mentioned_or('!'),
        intents = intents,
    ) as bot:
        await bot.start(token)
        await bot.load_cog('bot.cogs.my_cog')
asyncio.run(main())
