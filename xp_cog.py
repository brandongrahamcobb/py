""" xp_cog.py
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

import discord
from discord.ext import commands, tasks
import json
import random
import os

# Your XP table
xp_table = [1, 15, 34, 57, 92, 135, 372, 560, 840, 1242, 1144, 1573, 2144, 2800, 3640, 4700, 5893, 7360, 9144, 11120, 
            13477, 16268, 19320, 22880, 27008, 31477, 36600, 42444, 48720, 55813, 63800, 86784, 98208, 110932, 124432, 
            139372, 155865, 173280, 192400, 213345, 235372, 259392, 285532, 312928, 342624, 374760, 408336, 445544, 
            483532, 524160, 567772, 598886, 631704, 666321, 702836, 741351, 781976, 824828, 870028, 917625, 967995, 
            1021041, 1076994, 1136013, 1198266, 1263930, 1333194, 1406252, 1483314, 1564600, 1650340, 1740778, 1836173, 
            1936794, 2042930, 2154882, 2272970, 2397528, 2528912, 2667496, 2813674, 2967863, 3130502, 3302053, 3483005, 
            3673873, 3875201, 4087562, 4311559, 4547832, 4797053, 5059931, 5337215, 5629694, 5938202, 6263614, 6606860, 
            6968915, 7350811, 7753635, 8178534, 8626718, 9099462, 9598112, 10124088, 10678888, 11264090, 11881362, 
            12532461, 13219239, 13943653, 14707765, 15513750, 16363902, 17260644, 18206527, 19204245, 20256637, 21366700, 
            22537594, 23772654, 25075395, 26449526, 27898960, 29427822, 31040466, 32741483, 34535716, 36428273, 38424542, 
            40530206, 42751262, 45094030, 47565183, 50171755, 52921167, 55821246, 58880250, 62106888, 65510344, 69100311, 
            72887008, 76881216, 81094306, 85594273, 90225770, 95170142, 100385466, 105886589, 111689174, 117809740, 
            124265714, 131075474, 138258410, 145834970, 153826726, 162256430, 171148082, 180526997, 190419876, 200854885, 
            211861732, 223471711, 223471711, 248635353, 262260570, 276632449, 291791906, 307782102, 324648562, 342439302, 
            361204976, 380999008, 401877754, 423900654, 447130410, 471633156, 497478653, 524740482, 553496261, 583827855, 
            615821622, 649568646, 685165008, 722712050, 762316670, 804091623, 848155844, 894634784, 943660770, 995373379, 
            1049919840, 1107455447, 1168144006, 1232158297, 1299680571, 1370903066, 1446028554, 1525246918, 1608855764, 
            1697021059]

class XPCog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.user_data = self.load_data()
        self.announce_xp.start()

    def load_data(self):
        try:
            with open('../json/xp_data.json', 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            return {}

    def save_data(self):
        with open('../json/xp_data.json', 'w') as f:
            json.dump(self.user_data, f, indent=4)

    def get_level_from_xp(self, xp):
        """ Calculate the level based on XP. """
        level = 0
        while level < len(xp_table) and xp >= xp_table[level]:
            level += 1
        return level

    async def initialize_all_users(self):
        """ Initialize XP data for all current members in all guilds. """
        for guild in self.bot.guilds:
            for member in guild.members:
                if not member.bot:
                    user_id = str(member.id)
                    if user_id not in self.user_data:
                        self.user_data[user_id] = {'xp': 0, 'level': 0}
        self.save_data()

    @commands.Cog.listener()
    async def on_ready(self):
        print(f'Logged in as {self.bot.user}')
        await self.initialize_all_users()

    @commands.Cog.listener()
    async def on_member_join(self, member):
        user_id = str(member.id)
        if user_id not in self.user_data:
            self.user_data[user_id] = {'xp': 0, 'level': 0}
            self.save_data()
            print(f"Initialized XP data for new member: {member.name}")

    @commands.Cog.listener()
    async def on_member_remove(self, member):
        user_id = str(member.id)
        if user_id in self.user_data:
            del self.user_data[user_id]
            self.save_data()
            print(f"Reset XP data for member: {member.name}")

    @commands.Cog.listener()
    async def on_message(self, message):
        if message.author.bot:
            return

        user_id = str(message.author.id)
        if user_id not in self.user_data:
            self.user_data[user_id] = {'xp': 0, 'level': 0}

        current_xp = self.user_data[user_id]['xp']
        current_level = self.user_data[user_id]['level']

        # Determine the max possible XP increase based on current level
        if current_level < len(xp_table):
            max_increase = (xp_table[current_level] - current_xp) // 10
        else:
            max_increase = 50  # Cap the XP increase for very high levels

        # Add random XP amount
        xp_increase = random.randint(1, max_increase)
        self.user_data[user_id]['xp'] += xp_increase
        new_level = self.get_level_from_xp(self.user_data[user_id]['xp'])
        self.user_data[user_id]['level'] = new_level

        if new_level > current_level:
            await message.channel.send(f"{message.author.mention} has leveled up to level {new_level}!")

        self.save_data()

    @tasks.loop(minutes=5)
    async def announce_xp(self):
        channel = discord.utils.get(self.bot.get_all_channels(), name='bot-spam')
        if channel:
            sorted_users = sorted(self.user_data.items(), key=lambda x: x[1]['level'], reverse=True)
            leaderboard = "\n".join([f"<@{user_id}> - Level: {data['level']}, XP: {data['xp']}" for user_id, data in sorted_users[:10]])
            await channel.send(f"**XP Leaderboard:**\n{leaderboard}")

    @announce_xp.before_loop
    async def before_announce_xp(self):
        await self.bot.wait_until_ready()

    @commands.command(name='xp')
    async def xp(self, ctx, member: discord.Member = None):
        member = member or ctx.author
        user_id = str(member.id)
        if user_id not in self.user_data:
            self.user_data[user_id] = {'xp': 0, 'level': 0}

        xp = self.user_data[user_id]['xp']
        level = self.user_data[user_id]['level']

        await ctx.send(f'{member.mention} - Level: {level}, XP: {xp}')

    @commands.command(name='addxp')
    @commands.has_permissions(administrator=True)
    async def addxp(self, ctx, amount: int, member: discord.Member = None):
        member = member or ctx.author
        user_id = str(member.id)
        if user_id not in self.user_data:
            self.user_data[user_id] = {'xp': 0, 'level': 0}

        self.user_data[user_id]['xp'] += amount
        self.user_data[user_id]['level'] = self.get_level_from_xp(self.user_data[user_id]['xp'])
        self.save_data()

        await ctx.send(f'Added {amount} XP to {member.mention}. They are now level {self.user_data[user_id]["level"]}.')

    @commands.command(name='setlevel')
    @commands.has_permissions(administrator=True)
    async def setlevel(self, ctx, level: int, member: discord.Member = None):
        member = member or ctx.author
        user_id = str(member.id)
        if user_id not in self.user_data:
            self.user_data[user_id] = {'xp': 0, 'level': 0}

        self.user_data[user_id]['level'] = level
        self.user_data[user_id]['xp'] = xp_table[level - 1] if level > 0 else 0
        self.save_data()

        await ctx.send(f'Set {member.mention} to level {level}.')

async def setup(bot):
    await bot.add_cog(XPCog(bot))
