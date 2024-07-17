# cogs/message.py

import discord

from discord.ext import commands

class Message(commands.Cog):

    def __init__(self, bot):
        self.bot = bot
        self._last_member = None

    @commands.Cog.listener()
    @commands.is_owner()
    async def on_message(self, message: discord.Message) -> None:
        if message.author.bot:
            return

async def setup(bot):
    await bot.add_cog(Message(bot))
