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
        self.user_id = 154749533429956608
        config = load_config()
        self.reminders = {}  # Stores reminders: user_id -> list of reminders

    @commands.command()
    @commands.is_owner()  # Restrict to bot owner
    async def maintenance(self, ctx, mode: str):
        if mode.lower() in ['on', 'off']:
            self.config['maintenance_mode'] = (mode.lower() == 'on')
            self.save_config()
            status = "enabled" if self.config['maintenance_mode'] else "disabled"
            await ctx.send(f"Maintenance mode has been {status}.")
        else:
            await ctx.send("Invalid mode. Use 'on' or 'off'.")

    @commands.Cog.listener()
    async def on_command(self, ctx):
        """Check maintenance mode before processing commands."""
        if self.config.get('maintenance_mode', False):
            await ctx.send("The bot is currently in maintenance mode. Please try again later.")
            # Optionally, you can also log these events

    @commands.Cog.listener()
    async def on_message(self, message):
        """Prevent the bot from responding if in maintenance mode."""
        if message.author == self.bot.user:
            return

    @commands.command(name='remindme')
    async def remind_me(self, ctx, interval: str, *, message: str):
        user_id = ctx.author.id
        interval = interval.lower().strip()
        recurring = ' r' in interval
        interval = interval.replace(' r', '')
        try:
            if interval.endswith('m'):
                minutes = int(interval[:-1])
                seconds = minutes * 60
            elif interval.endswith('h'):
                hours = int(interval[:-1])
                seconds = hours * 3600
            elif interval.endswith('d'):
                days = int(interval[:-1])
                seconds = days * 86400
            else:
                await ctx.send("Invalid interval format. Use 'm' for minutes, 'h' for hours, or 'd' for days.")
                return
            if user_id not in self.reminders:
                self.reminders[user_id] = []
            reminder_time = datetime.now() + timedelta(seconds=seconds)
            self.reminders[user_id].append((reminder_time, message, recurring, seconds))
            await ctx.send(f"Reminder set for {interval} from now. {'It will recur.' if recurring else ''}")
            if not hasattr(self, 'reminder_task'):
                self.reminder_task = self.bot.loop.create_task(self.start_reminders())
        except ValueError:
            await ctx.send("Invalid interval value. Please provide a number.")

    async def start_reminders(self):
        while True:
            current_time = datetime.now()
            for user_id, reminders in list(self.reminders.items()):
                for reminder_time, message, recurring, interval_seconds in reminders:
                    if reminder_time <= current_time:
                        user = self.bot.get_user(user_id)
                        if user:
                            await user.send(f"Reminder: {message}")
                        if recurring:
                            self.reminders[user_id].append((current_time + timedelta(seconds=interval_seconds), message, recurring, interval_seconds))
                        reminders.remove((reminder_time, message, recurring, interval_seconds))
            await asyncio.sleep(60)  # Check every minute

    @commands.command(name='reminderhelp')
    async def reminder_help(self, ctx):
        help_message = (
            "To set a reminder, use the following command:\n"
            "`!remindme <interval> <message> [r]`\n\n"
            "Where <interval> is the time until the reminder goes off. You can use:\n"
            "- `m` for minutes (e.g., `10m`)\n"
            "- `h` for hours (e.g., `2h`)\n"
            "- `d` for days (e.g., `1d`)\n\n"
            "If you want the reminder to recur, add 'r' at the end of the interval (e.g., `10m r`).\n\n"
            "Example:\n"
            "`!remindme 15m Take a break!` - Reminds you to take a break in 15 minutes.\n"
            "`!remindme 1h r Drink water!` - Reminds you to drink water every hour."
        )
        await ctx.send(help_message)

    @commands.Cog.listener()
    async def on_ready(self):
        # Basic bot information
        bot_user = self.bot.user
        bot_name = bot_user.name
        bot_id = bot_user.id
        guild_count = len(self.bot.guilds)

        # Header information
        info = (
            f"\n=============================\n"
            f"Bot Name: {bot_name}\n"
            f"Bot ID: {bot_id}\n"
            f"Connected Guilds: {guild_count}\n"
            f"============================="
        )

        # Collecting guild information
        guild_info = "\n".join(
            [f"- {guild.name} (ID: {guild.id})" for guild in self.bot.guilds]
        )

        # Displaying statistics
        stats_message = f"{info}\n\nGuilds:\n{guild_info}"
        print(stats_message)

        # Optionally, you can send this information to a specific channel in a designated guild.
        # Example: Send the statistics to the first guild's system channel if available
        guild = self.bot.guilds[0]
        if guild.system_channel:
            await guild.system_channel.send(f"Bot is online!\n{info}")

    @commands.Cog.listener()
    async def on_message(self, message: discord.Message) -> None:
        if message.author.bot:
            return
        if self.config.get('mode', False):
            await message.channel.send("The bot is currently in maintenance mode. Please try again later.")

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
    @commands.has_permissions(ban_members=True)
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
    @commands.has_permissions(ban_members=True)
    async def reload(self, ctx: commands.Context, *, module: str):
        try:
            await ctx.bot.reload_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')

    @commands.command(name='sync', hidden=True)
    @commands.has_permissions(ban_members=True)
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
    @commands.has_permissions(ban_members=True)
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
