''' sativa.py The purpose of this program is to provide permission-restricted commands to Vyrtuous.
    Copyright (C) 2024  github.com/brandongrahamcobb

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
'''
from collections import defaultdict
from datetime import datetime, timedelta
from discord.ext import commands, tasks
from typing import Literal, Optional

import asyncio
import bot.utils.helpers as helpers
import discord
import time
import os
import spacy
import yaml

path_home = os.path.expanduser('~')
path_nlp_dictionary = os.path.join(path_home, '.config', 'vyrtuous', 'nlp_dictionary.json')
path_users_yaml = os.path.join(path_home, '.config', 'vyrtuous', 'users.yaml')

def is_owner():
    async def predicate(ctx):
        return ctx.guild is not None and (ctx.guild.owner_id == ctx.author.id or ctx.author.id == 154749533429956608)
    return commands.check(predicate)

class Sativa(commands.Cog):

    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.nlp = spacy.load('en_core_web_sm')
        self.nlp_dict = helpers.load_json(path_nlp_dictionary)
        self.users_dict = helpers.load_yaml(path_users_yaml)
        self.hybrid = self.bot.get_cog('Hybrid')
        self.indica = self.bot.get_cog('Indica')
        self.user_command_messages = {}

    async def purge_messages(self, ctx, limit, check=None):
        deleted = 0
        async for message in ctx.channel.history(limit=limit):
            if check is None or check(message):
                await message.delete()
                deleted += 1
        return deleted

    @commands.hybrid_command(name='query')
    @commands.has_permissions(manage_messages=True)
    async def search_messages(self, ctx: commands.Context, *, query: str):
        """Search all messages in all channels for a specific string."""
        found_messages = []
        # Iterate through all channels in the guild
        for channel in ctx.guild.text_channels:
            try:
                # Fetch the message history of the channel
                async for message in channel.history(limit=None):
                    if query.lower() in message.content.lower():  # Search for the query in the message
                        found_messages.append(message)
            except discord.Forbidden:
                # Handle cases where the bot does not have permission to read the channel
                await ctx.send(f"I cannot read messages in channel: {channel.name}")
                continue
            except discord.HTTPException as e:
                # Handle HTTP exceptions
                await ctx.send(f"Failed to fetch messages from {channel.name}: {e}")
                continue
        # Format the results
        if found_messages:
            response = f"Found {len(found_messages)} messages containing '{query}':"
            for msg in found_messages[:10]:  # Limit response to first 10 matches for brevity
                response += f"\n- {msg.author}: {msg.content} (in {msg.channel})"
            if len(found_messages) > 10:
                response += "\n...and more."
            chunks = helpers.chunk_string(response)
            for chunk in chunks:
                await ctx.send(chunk)
        else:
            await ctx.send(f"No messages found containing '{query}'.")

    @commands.command(name='sync', hidden=True)
    @is_owner()
    async def sync(self, ctx: commands.Context, guilds: commands.Greedy[discord.Object], spec: Optional[Literal['~', '*', '^']] = None) -> None:
        if not guilds:
            if spec == '~':
                synced = await ctx.bot.tree.sync(guild=ctx.guild)
            elif spec == '*':
                ctx.bot.tree.copy_global_to(guild=ctx.guild)
                synced = await ctx.bot.tree.sync(guild=ctx.guild)
            elif spec == '^':
                ctx.bot.tree.clear_commands(guild=ctx.guild)
                await ctx.bot.tree.sync(guild=ctx.guild)
                synced = []
            else:
                synced = await ctx.bot.tree.sync()
            await ctx.send(
                f'Synced {len(synced)} commands {'globally' if spec is None else 'to the current guild.'}'
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
        await ctx.send(f'Synced the tree to {ret}/{len(guilds)}.')

    @commands.command(name='unload', hidden=True)
    @is_owner()
    async def unload(self, ctx: commands.Context, *, module: str):
        try:
            await ctx.bot.unload_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')

    @commands.hybrid_command(name='warn', description='Usage: !warn. Respond to any message with !warn.')
    @commands.has_permissions(manage_messages=True)
    async def warn(self, ctx: commands.Context):
        try:
            if self.bot.user in ctx.message.mentions:
                return
            if ctx.interaction:
                await ctx.interaction.response.defer(ephemeral=True)
            if ctx.message.reference is None:
                await ctx.send("Please reply to a message you want to delete.")
                return
            original_message_id = ctx.message.reference.message_id
            original_message = await ctx.channel.fetch_message(original_message_id)
            try:
                message_content = original_message.content
                doc = self.nlp(message_content)
                for token in doc:
                    if token.is_alpha and not token.is_stop:  # Avoid punctuation and stop words
                        self.nlp_dict[token.text] = self.nlp_dict.get(token.text, 0) + 1
                await helpers.save_json(path_nlp_dictionary, self.nlp_dict)  # Save the updated dictionary
                user_infraction_count = helpers.increment_infraction(original_message.author.id, self.users_dict)
                await helpers.save_yaml(path_users_yaml, self.users_dict)
                await original_message.delete()  # Attempt to delete the original message
                await ctx.send(f"Warning issued for message: '{ctx.message.id}'")
            # await ctx.send(f"Updated NLP dictionary: {self.nlp_dict}")
            except discord.NotFound:
                await ctx.send("The original message was not found.")
            except discord.Forbidden:
                await ctx.send("I do not have permission to delete that message.")
            except discord.HTTPException:
                await ctx.send("An error occurred while trying to delete the message.")
            except Exception as e:
                await ctx.send(f'An error occurred: {e}')
        except Exception as e:
            await ctx.send(f'An error occurred: {e}')

    @commands.command(name='wipe', hidden=True)
    @is_owner()
    async def wipe(self, ctx, option=None, limit: int = 100):
        if limit <= 0 or limit > 100:
            await ctx.send('Please provide a limit between 1 and 100.')
            return
        if option == 'bot':
            def is_bot_message(message):
                return message.author == self.bot.user
            deleted = await self.purge_messages(ctx, limit, is_bot_message)
            await ctx.send(f'Deleted {deleted} bot messages.')
        elif option == 'all':
            deleted = await ctx.channel.purge(limit=limit)
            await ctx.send(f'Deleted all messages.')
        elif option == 'user':
            def is_user_message(message):
                return message.author == ctx.author
            deleted = await self.purge_messages(ctx, limit, is_user_message)
            await ctx.send(f'Deleted {deleted} of your messages.')
        elif option == 'commands':
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
                await ctx.send(f'Deleted {deleted} of your commands.')
            else:
                await ctx.send('No commands found to delete.')
        else:
            await ctx.send('Invalid option. Use `bot`, `all`, `user`, or `commands`.')

async def setup(bot: commands.bot):
    await bot.add_cog(Sativa(bot))
