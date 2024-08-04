""" game_cog.py
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
from discord.ext import commands
import random
import json
import os

class GameCog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.users = {}
        self.load_users()
        self.xp_table = [1, 15, 34, 57, 92, 135, 372, 560, 840, 1242, 1144, 1573, 2144, 2800, 3640, 4700, 5893, 7360, 9144, 11120, 
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

    def load_users(self):
        if os.path.exists("users.json"):
            with open("users.json", "r") as f:
                self.users = json.load(f)
    
    def save_users(self):
        with open("users.json", "w") as f:
            json.dump(self.users, f, indent=4)

    def get_xp_for_level(self, level):
        return sum(self.xp_table[:level])

    def get_xp_per_interaction(self, current_level):
        xp_for_next_level = self.xp_table[current_level]
        daily_xp_required = xp_for_next_level / (365 / 200)
        xp_per_interaction = daily_xp_required / 200
        return xp_per_interaction

    def distribute_xp(self, user_id):
        if user_id not in self.users:
            self.users[user_id] = {
                "xp": 0,
                "level": 1
            }
        current_level = self.users[user_id]["level"]
        xp_per_interaction = self.get_xp_per_interaction(current_level)
        xp_gain = random.uniform(0.8 * xp_per_interaction, 1.2 * xp_per_interaction)
        self.users[user_id]["xp"] += xp_gain

        if self.users[user_id]["xp"] >= self.get_xp_for_level(current_level + 1):
            self.users[user_id]["level"] += 1

        self.save_users()

    @commands.Cog.listener()
    async def on_message(self, message):
        if message.author.bot:
            return

        self.distribute_xp(str(message.author.id))

        if self.users[str(message.author.id)]["xp"] >= self.get_xp_for_level(self.users[str(message.author.id)]["level"] + 1):
            await message.channel.send(f"Congratulations {message.author.mention}, you reached level {self.users[str(message.author.id)]['level']}!")

    @commands.command()
    async def level(self, ctx, member: discord.Member = None):
        member = member or ctx.author
        user_id = str(member.id)
        self.distribute_xp(user_id)  # Distribute XP before displaying level and XP

        if user_id in self.users:
            level = self.users[user_id]["level"]
            xp = self.users[user_id]["xp"]
            xp_needed_for_next_level = self.get_xp_for_level(level + 1) - xp
            await ctx.send(f"{member.mention} is at level {level} with {xp:.2f} XP.\n"
                           f"You need {xp_needed_for_next_level:.2f} XP to reach level {level + 1}.")
        else:
            await ctx.send(f"{member.mention} has not interacted with the bot yet.")

    @commands.command()
    async def leaderboard(self, ctx):
        # Distribute XP to the user invoking the command
        self.distribute_xp(str(ctx.author.id))

        # Sort users by level, then by XP
        sorted_users = sorted(self.users.items(), key=lambda item: (item[1]["level"], item[1]["xp"]), reverse=True)
        
        # Create leaderboard message
        leaderboard_message = "Leaderboard:\n"
        for i, (user_id, data) in enumerate(sorted_users[:10], start=1):
            member = ctx.guild.get_member(int(user_id))
            member_name = member.display_name if member else "Unknown User"
            level = data["level"]
            xp = data["xp"]
            leaderboard_message += f"{i}. {member_name} - Level {level}, {xp:.2f} XP\n"

        await ctx.send(leaderboard_message)

async def setup(bot):
    await bot.add_cog(GameCog(bot))
