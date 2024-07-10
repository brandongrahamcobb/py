# Listener.py
import discord
import os
from discord.ext import commands
from subprocess import run

class Listener(commands.Cog):

    def __init__(self, bot):
        self.bot = bot
        self_last_member = None

    def version():
        if not os.path.exists('../txt/version.txt'):
             with open('../txt/version.txt', 'w') as f:
                 f.write('1.0.0')
        with open('../txt/version.txt', 'r') as f:
             version = f.read().strip()
        major, minor, patch = map(int, version.split('.'))
        patch += 1
        if patch >= 10:
            patch = 0
            minor += 1
        if minor >= 10:
            minor = 0
            major += 1
        version = f"{major}.{minor}.{patch}"
        with open('../txt/version.txt', 'w') as f:
            f.write(version)
        return version

    @commands.Cog.listener()
    @commands.is_owner()
    async def on_message(self, message: discord.Message) -> None:
        if message.author.bot:
            return

    @commands.Cog.listener()
    async def on_ready(self):
        version = Listener.version()
        print('Logged on as', self.bot.user)
        print(f'Bot ID: {self.bot.user.id}')
        print(f'Bot Version: {version}')
        print(f'Connected to {len(self.bot.guilds)} guilds:')
        for guild in self.bot.guilds:
            print(f' - {guild.name} (id: {guild.id})')

async def setup(bot):
    await bot.add_cog(Listener(bot))
