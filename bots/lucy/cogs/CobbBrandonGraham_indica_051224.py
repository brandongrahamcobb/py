import json
import discord
from discord.ext import commands
import bot.utils.helpers as helpers
import os

path_home = os.path.expanduser('~')
path_users_yaml = os.path.join(path_home, '.config', 'vyrtuous', 'users.yaml')

class Indica(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.users = {}  # Assuming this is initialized with your user data structure

    @commands.Cog.listener()
    async def on_message(self, message):
        if message.author == self.bot.user:
            return  # Ignore messages from the bot itself

        if self.bot.user.mentioned_in(message) or isinstance(message.channel, discord.DMChannel):
            moderation_input = f'{message.content} in channel: {message.channel.name}'
            structured_sys_input = self.get_sys_input(message.id)
            try:
                # Assuming deprecated_create_moderation is an asynchronous generator
                async for moderation_response in helpers.deprecated_create_moderation(
                    input_text=moderation_input,
                    sys_input=structured_sys_input,
                    conversation_id='1'
                ):
                    # Log the raw response for debugging
                    print("Moderation Response Raw Output:", repr(moderation_response))

                    # Validate the response before parsing
                    if not moderation_response or not moderation_response.strip():
                        print("Received an empty or null response from moderation.")
                        continue  # Skip this iteration and continue with the next

                    # Attempt to parse the moderation response
                    try:
                        moderation_dict = json.loads(moderation_response)
                    except json.JSONDecodeError as e:
                        print(f"JSON decoding failed: {e}")
                        print(f"Response content was: '{moderation_response}'")
                        # Handle or log the error as needed
                        continue  # Skip this iteration and continue with the next

                    # Perform further processing with the valid JSON
                    flagged = moderation_dict.get('results', [{}])[0].get('flagged', False)
                    await self.handle_moderation_result(flagged, message)
            except Exception as e:
                print(f"Error during moderation: {e}")

    def get_sys_input(self, message_id):
        # Provide clear instructions for both input and expected output
        return f"""
You are a moderation assistant. Please analyze the following message and respond in this structured JSON format:
{{
    "id": "{message_id}",
    "model": "gpt-4o-mini",
    "results": [
        {{
            "flagged": false,
            "categories": {{
                "sexual": false,
                "sexual/minors": false,
                "harassment": false,
                "harassment/threatening": false,
                "violence": false,
                "violence/graphic": false,
                "self-harm": false,
                "self-harm/intent": false,
                "self-harm/instructions": false,
                "self-harm/ideation": false,
                "hate": false,
                "hate/threatening": false,
                "extremism": false,
                "extremism/violence": false
            }},
            "category_scores": {{
                "sexual": 0,
                "sexual/minors": 0,
                "harassment": 0,
                "harassment/threatening": 0,
                "violence": 0,
                "violence/graphic": 0,
                "self-harm": 0,
                "self-harm/intent": 0,
                "self-harm/instructions": 0,
                "self-harm/ideation": 0,
                "hate": 0,
                "hate/threatening": 0,
                "extremism": 0,
                "extremism/violence": 0
            }}
        }}
    ]
}}
"""
    async def handle_moderation_result(self, flagged, message):
        if flagged:
            await message.delete()
            count = helpers.increment_infraction(str(message.author.id), self.users)
            await helpers.save_yaml(path_users_yaml, self.users)
            embed = helpers.create_embed(
                "Message Deleted",
                f"Your message was deleted due to inappropriate content. Infractions: {count}"
            )
            await message.author.send(embed=embed)
        else:
            # This may also call for a response
            response = await helpers.deprecated_create_completion(f'{message.content}', self.get_sys_input(message.id))
            responses = helpers.chunk_string(text=response)
            await message.reply(f"@{message.author.name}, {responses[0]}")
            for response in responses[1:]:
                await message.reply(response)

    @commands.Cog.listener()
    async def on_message_edit(self, before, after):
        if before.content != after.content:
            ctx = await self.bot.get_context(after)
            if ctx.command:
                await self.bot.invoke(ctx)

    @commands.Cog.listener()
    async def on_command_error(self, ctx: commands.Context, error):
        # Handle command errors properly
        if isinstance(error, commands.CommandNotFound):
            await ctx.send("Unknown command. Use `!help` to see available commands.")
        elif isinstance(error, commands.MissingRequiredArgument):
            await ctx.send("Missing arguments. Use `!help` to see how to use the command.")
        elif isinstance(error, commands.BadArgument):
            await ctx.send("Invalid argument type. Please check your inputs.")
        else:
            await ctx.send(f"An error occurred: {str(error)}")

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

async def setup(bot: commands.Bot):
    await bot.add_cog(Indica(bot))
