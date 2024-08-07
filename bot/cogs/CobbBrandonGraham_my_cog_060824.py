""" my_cog.py
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
import discord
import traceback

class MyCog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot

    @commands.command(hidden='true')
    @commands.is_owner()
    async def test(self, ctx, *args):
        try:
            await ctx.send(args)
        except Exception as e:
            await ctx.send(f'Error fetching data: {e}')

async def setup(bot: commands.Bot):
    await bot.add_cog(MyCog(bot))
