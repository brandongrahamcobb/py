# Super.py

from collections import Counter
from discord.ext import commands
from PIL import Image, ImageDraw

import discord
import io
import os
import requests

class Super(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self._last_member = None

    @commands.command(description='')
    async def color(self, ctx: commands.Context, *args):
        r = int(args[0])
        g = int(args[1])
        b = int(args[2])
        guildroles = await ctx.guild.fetch_roles()
        position = len(guildroles) - 1
        for arg in ctx.author.roles:
            if arg.name.isnumeric():
                await ctx.author.remove_roles(arg)
        for arg in guildroles:
            if arg.name.lower() == f'{r}{g}{b}':
                await ctx.author.add_roles(arg)
                await arg.edit(position=position)
                await ctx.send(f'I successfully changed your role color to {r}, {g}, {b}')
                return
        newrole = await ctx.guild.create_role(name=f'{r}{g}{b}', color=discord.Color.from_rgb(r, g, b), reason='new color')
        await newrole.edit(position=position)
        await ctx.author.add_roles(newrole)
        await ctx.send(f'I successfully changed your role color to {r}, {g}, {b}')

async def setup(bot):
    await bot.add_cog(Super(bot))

