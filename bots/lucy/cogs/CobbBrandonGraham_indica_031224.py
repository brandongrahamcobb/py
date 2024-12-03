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
from discord.ext import commands, tasks
import bot.utils.helpers as helpers

import asyncio
import discord

class Indica(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.hybrid = self.bot.get_cog('Hybrid')
        self.sativa = self.bot.get_cog('Sativa')
        self.sys_input = """
        """
        self.channel_id = 1305608084017905734
#        self.post_message.start()

    def cog_unload(self):
        self.post_message.stop()

    @tasks.loop(hours=5)  # Set the loop to run every 5 minutes
    async def post_message(self):
        channel = self.bot.get_channel(self.channel_id)
        if channel is not None:
            async for response in helpers.deprecated_create_completion('Inspire me to be hungry!', self.sys_input, self.channel_id):
                responses = helpers.chunk_string(text=response)
                for response in responses:
                    await channel.send(response)

    @post_message.before_loop
    async def before_post_message(self):
        await self.bot.wait_until_ready()  # Wait until the bot is ready
#            Do not use idioms which use animals as part of the idiom.
#            Never suggest animal products.
#            Forified foods are not a good source of B12.
#            Advocate for B12 from seaweed.
#            I want you to never accept vegetarianism, it must be fully vegan.
#            You are an AI you are uneqipped to be giving out information about drugs
#            I want you to verify whether the user can answer these following questions:
#            Name 5 reasons that eating animal products is unsustainable? 
#            Why is sterilization bad for a sustainable microbiome?
#            Why are antibiotics bad for a sustainable microbiome
#            Why is an animal product based diet unsustainable 
#            Why is a high protein diet unsustainable?
#            Why is a high sodium diet unsustainable
#            Why is a diet high in processed sugar unsustainable 
#            A plant base diet can help prevent and reverse disease. T of F 
#            Whole foods arent necessary for a sustainable diet. T or F
#            Why isnt supplementation a sustainable lifestyle choice?
#            Why is a high fat diet unsustainable?
#            Why is it unsustainable to eat a diet filled with processed foods
#            Why is inflammation important. Name five unsustainable inflammatory habits
#            Why is the gut-brain connection so important to the immune system (vagus nerve bundle, nervous system)
#            Name 3 gut connection axis\'s (ex. gut-brain) 
#            Name 7 myths about supplemental requirements on a plant based diet
#            Describe ways that the nervous system is affected by diet choices
#            Why is fiber important for diet
#            Why is the gylcemic index important
#            What are phytochemicals, and why are they important to health
#            Name 3 reasons a low carb diet is unsustainable . Alternately, describe why low carb diets are unsustainable 
#            What is vitamin d? 
#            Can type 2 diabetes be reversed? 
#            Can pandemics be avoided? 
#            What is precision fermentation and why is it important?
#            *..is hibernation sustainable?
#            Why  is regenerative farming unsustainable?
#            Explain why sea vegetables are so valuable to a nutritious diet
#            Is added salt dangerous?
#            Why is acidity important to diet?
#            Why does wealth inequality matter? 
#            why are t cells so important to immunology?
#            Why is sustainable hygiene important? 
#            Why is the mesolimbic system important in addiction?
#            Why arent supplements sustainable?
#            Why are phytonutrients important?
#            Why is a calorie, calorie out a short sighted approach to nutrition

#    async def check_message_for_misinformation(self, message):
#        await message.channel.send(f"⚠️ Misinformation alert: {response}")

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
        if message.author == self.bot.user or '!' in message.content[0]:
            return
        if self.bot.user.mentioned_in(message) or isinstance(message.channel, discord.DMChannel):
            if int(message.guild.id) != int(self.config['testing_guild_id']):
                conversation_id = message.channel.id + message.author.id
                async for response in helpers.deprecated_create_completion(f'{message.content}', self.sys_input, conversation_id):
                    responses = helpers.chunk_string(text=response)
                    for response in responses:
                        await message.channel.send(f"@{message.author.name}, {response}")
        if int(message.guild.id) == int(self.config['testing_guild_id']):
            conversation_id = None
            async for response in helpers.deprecated_create_completion(message.content, sys_input="Verify whether the prompt is about {message.channel.name}. If it isn't, respond. Otherwise, respond with an empty string", conversation_id=None):
                if response != '':
                    async with message.channel.typing():
                        await message.channel.send(response)
                        await message.delete()
#        response = await self.check_for_misinformation(message.content)
 #       if response:
  #          await message.channel.send(f"@{message.author.name}, {response}")
   #         async for response in helpers.deprecated_create_completion(prompt, self.sys_input, "response"):
    #            if response:
     #               async with message.channel.typing():
      #                  await asyncio.sleep(1)
       #                 await message.channel.send(response)
#        await self.bot.process_commands(message)

        
   #     async with message.channel.typing():
        #    async for response in helpers.create_completion(f'{message.content}', message.author.id):
        # Spoiler content-warning material
#        if message.author.id in self.users:
#            spoiler_content = f"||{message.content}||"
#            await message.delete()
#            await message.channel.send(f"{message.author.mention} said: {spoiler_content}")

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
        user = await self.bot.fetch_user(self.config['owner_id'])
        await user.send(stats_message)

async def setup(bot: commands.Bot):
    await bot.add_cog(Indica(bot))
