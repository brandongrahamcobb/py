# cogs/kick.py

import discord

from discord.ext import commands

class Kick(commands.Cog):

    def __init__(self, bot):
        self.bot = bot
        self._last_member = None

    @commands.command()
    @commands.is_owner()
    async def kick(self, ctx: commands.Context, member: discord.Member, *, reason: str):
        await member.kick(reason=reason)
        await ctx.send(f'Kicked {member.name} for: {reason}')

async def setup(bot):
    await bot.add_cog(Kick(bot))
