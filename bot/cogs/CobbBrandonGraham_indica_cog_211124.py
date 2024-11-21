''' ListenerCog.py The purpose of the program is to contain a declaration of a few self variables and all the listeners for a discord bot.
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
from discord.ext import commands
import bot.utils.helpers as helpers

import discord

class Indica(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.hybrid = self.bot.get_cog('Hybrid')
        self.sativa = self.bot.get_cog('Sativa')

    @commands.Cog.listener()
    async def on_message_edit(self, before, after):
        if before.content != after.content:
            ctx = await self.bot.get_context(after)
            if ctx.command:
                await self.bot.invoke(ctx)

    @commands.Cog.listener()
    async def on_member_join(self, member = discord.Member):
        await member.send("DM me to chat!")

    @commands.Cog.listener()
    async def on_message(self, message):
        # Ignore messages from self.
        if message.author == self.bot.user:
            return
        # Intepret messages for AI
        if self.bot.user.mentioned_in(message) or isinstance(message.channel, discord.DMChannel):
            conversation_id = message.author.id
            async with message.channel.typing():
                async for response in helpers.create_completion(f'{message.content}', message.author.id):
                    await message.channel.send(response)
        # Spoiler content-warning material
        if message.author.id in self.users:
            spoiler_content = f"||{message.content}||"
            await message.delete()
            await message.channel.send(f"{message.author.mention} said: {spoiler_content}")

    @commands.Cog.listener()
    async def on_ready(self):
        bot_user = self.bot.user
        bot_name = bot_user.name
        bot_id = bot_user.id
        guild_count = len(self.bot.guilds)
        info = (
            f'\n=============================\n'
            f'bot Name: {bot_name}\n'
            f'bot ID: {bot_id}\n'
            f'Connected Guilds: {guild_count}\n'
            f'============================='
        )
        guild_info = '\n'.join(
            [f'- {guild.name} (ID: {guild.id})' for guild in self.bot.guilds]
        )
        stats_message = f'{info}\n\nGuilds:\n{guild_info}'
        print(stats_message)
        user = await self.bot.fetch_user(154749533429956608)
        await user.send(stats_message)

async def setup(bot: commands.Bot):
    await bot.add_cog(Indica(bot))
