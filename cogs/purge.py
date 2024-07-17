# cogs/Purge.py

import discord

from discord.ext import commands

class Purge(commands.Cog):

    def __init__(self, bot):
        self.bot = bot
        self._last_member = None

    @commands.command()
    @commands.is_owner()
    async def purge(self, ctx: commands.Context):
        def not_pinned(msg):
            return not msg.pinned
        purged = await ctx.channel.purge(limit=100, check=not_pinned)
        await ctx.send(f"Successfully removed {len(purged)} non-pinned messages!")

async def setup(bot):
    await bot.add_cog(Purge(bot))
