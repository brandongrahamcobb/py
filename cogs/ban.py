# cogs/ban.py

import discord

from discord.ext import commands

class Ban(commands.Cog):

    def __init__(self, bot):
        self.bot = bot
        self._last_member = None

    @commands.command()
    @commands.is_owner()
    async def ban(self, ctx: commands.Context, member: discord.Member, reason: str):
        await member.ban(reason=reason)
        await ctx.send(f'Banned {member.name} for: {reason}')

async def setup(bot):
    await bot.add_cog(Ban(bot))
