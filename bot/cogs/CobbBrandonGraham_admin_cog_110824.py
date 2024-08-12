""" admin_cog.py
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

from discord.ext import commands, tasks
from bot.utils.helpers import load_config

import asyncio
import discord
import json
import io
import os
import pubchempy as pcp
import re
import traceback

from datetime import datetime, timedelta
from typing import Literal, Optional

class AdminCog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        config = load_config()
        self.testing_guild_id = self.bot.testing_guild_id
        self.user_command_messages = {}

    @commands.Cog.listener()
    async def on_message(self, message):
        if message.author == self.bot.user:
            return
        if message.author != self.bot.user and message.content.startswith('!'):
            if message.author.id not in self.user_command_messages:
                self.user_command_messages[message.author.id] = []
            self.user_command_messages[message.author.id].append(message.id)

    @commands.Cog.listener()
    async def on_ready(self):
        bot_user = self.bot.user
        bot_name = bot_user.name
        bot_id = bot_user.id
        guild_count = len(self.bot.guilds)
        info = (
            f"\n=============================\n"
            f"Bot Name: {bot_name}\n"
            f"Bot ID: {bot_id}\n"
            f"Connected Guilds: {guild_count}\n"
            f"============================="
        )
        guild_info = "\n".join(
            [f"- {guild.name} (ID: {guild.id})" for guild in self.bot.guilds]
        )
        stats_message = f"{info}\n\nGuilds:\n{guild_info}"
        print(stats_message)

    @commands.Cog.listener()
    async def on_message(self, message: discord.Message) -> None:
        if message.author.bot:
            return

    @commands.Cog.listener()
    async def on_message_edit(self, before, after):
        if before.author == self.bot.user:
            return

    @commands.command(name='botinfo')
    @commands.is_owner()  # Ensures only the bot owner can use this command
    async def botinfo(self, ctx, member: discord.Member, prompt: str):
        """DMs the specified user with the bot's information."""
        bot_info = (
            "Here is some information about the bot:\n"
            f"**Name:** {self.bot.user.name}\n"
            f"**ID:** {self.bot.user.id}\n"
            f"**Description:** This is a bot created by Brandon Graham Cobb for various functionalities.\n"
            f"**Invite Link:** [Invite Bot](https://discord.com/oauth2/authorize?client_id=302202228016414721)\n"
            f"**Support Server:** [discord.gg](https://discord.gg/dNfUn8MeYN)"
            f"**Source:** [GitHub](https://github.com/brandongrahamcobb/py)"
        )
        try:
            await member.send(bot_info)
            await ctx.send(f"Sent bot information to {member.mention}.")
        except discord.Forbidden:
            await ctx.send(f"Failed to send DM to {member.mention}. They might have DMs disabled.")

    @commands.command(description='')
    @commands.has_permissions(manage_roles=True)
    async def colorize(self, ctx: commands.Context, *args):
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

    @commands.command()
    async def dmpurge(self, ctx: commands.Context):
        async for message in ctx.author.history(limit=9999999):
            if message.author.id == self.bot.user.id:
                await message.delete()

    @commands.command(name='load', hidden=True)
    @commands.is_owner()
    async def load(self, ctx: commands.Context, *, module: str):
        try:
            await ctx.bot.load_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')

    async def purge_messages(self, ctx, limit, check=None):
        deleted = 0
        async for message in ctx.channel.history(limit=limit):
            if check is None or check(message):
                await message.delete()
                deleted += 1
        return deleted

    @commands.hybrid_command()
    @commands.is_owner()
    async def reload(self, ctx: commands.Context, *, module: str):
        try:
            await ctx.bot.reload_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')

    @commands.command(name='sync', hidden=True)
    @commands.is_owner()
    async def sync(self, ctx: commands.Context, guilds: commands.Greedy[discord.Object], spec: Optional[Literal["~", "*", "^"]] = None) -> None:
        if not guilds:
            if spec == "~":
                synced = await ctx.bot.tree.sync(guild=ctx.guild)
            elif spec == "*":
                ctx.bot.tree.copy_global_to(guild=ctx.guild)
                synced = await ctx.bot.tree.sync(guild=ctx.guild)
            elif spec == "^":
                ctx.bot.tree.clear_commands(guild=ctx.guild)
                await ctx.bot.tree.sync(guild=ctx.guild)
                synced = []
            else:
                synced = await ctx.bot.tree.sync()
            await ctx.send(
                f"Synced {len(synced)} commands {'globally' if spec is None else 'to the current guild.'}"
            )
            return
        ret = 0
        for guild in guilds:
            try:
                await ctx.bot.tree.sync(guild=guild)
            except discord.HTTPException:
                pass
            else:
                ret += 1
        await ctx.send(f"Synced the tree to {ret}/{len(guilds)}.")

    @commands.command(name='unload', hidden=True)
    @commands.is_owner()
    async def unload(self, ctx: commands.Context, *, module: str):
        try:
            await ctx.bot.unload_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')

    @commands.command()
    @commands.has_permissions(manage_messages=True)
    async def wipe(self, ctx, option=None, limit: int = 100):
        if limit <= 0 or limit > 100:
            await ctx.send("Please provide a limit between 1 and 100.")
            return
        if option == "bot":
            def is_bot_message(message):
                return message.author == self.bot.user
            deleted = await self.purge_messages(ctx, limit, is_bot_message)
            await ctx.send(f"Deleted {deleted} bot messages.")
        elif option == "all":
            deleted = await ctx.channel.purge(limit=limit)
          #  deleted = await self.purge_messages(ctx, limit)
            await ctx.send(f"Deleted {deleted} messages.")
        elif option == "user":
            def is_user_message(message):
                return message.author == ctx.author
            deleted = await self.purge_messages(ctx, limit, is_user_message)
            await ctx.send(f"Deleted {deleted} of your messages.")
        elif option == "commands":
            if ctx.author.id in self.user_command_messages:
                message_ids = self.user_command_messages[ctx.author.id]
                deleted = 0
                for message_id in message_ids:
                    try:
                        message = await ctx.channel.fetch_message(message_id)
                        await message.delete()
                        deleted += 1
                    except discord.NotFound:
                        continue
                del self.user_command_messages[ctx.author.id]
                await ctx.send(f"Deleted {deleted} of your commands.")
            else:
                await ctx.send("No commands found to delete.")
        else:
            await ctx.send("Invalid option. Use `bot`, `all`, `user`, or `commands`.")

async def setup(bot: commands.Bot):
    await bot.add_cog(AdminCog(bot))
