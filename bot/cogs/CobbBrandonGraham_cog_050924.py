''' my_cog.py
    Copyright (C) 2024 github.com/<your-github>

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
from gradio_client import Client
from openai import AsyncOpenAI
from random import randint
import openai
import os

import discord

def is_owner():
    async def predicate(ctx):
       return ctx.guild is not None and (ctx.guild.owner_id == ctx.author.id or ctx.author.id == 154749533429956608)
    return commands.check(predicate)

class MyCog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.config = self.bot.config
        self.user_context = {}  # Dictionary to hold conversation history for each user

    def add_to_user_context(self, user_id, message):
        context = self.get_user_context(user_id)
        context.append(message)
        if len(context) > 10:  # Limit context length
            context.pop(0)

    def get_user_context(self, user_id):
        return self.user_context.setdefault(user_id, [])

    async def chatgpt(self, message: str):
        client = AsyncOpenAI(api_key=self.config['api_keys']['api_key_1'])
        try:
            stream = await client.chat.completions.create(
                model='gpt-4o-mini',
                messages=[
                    {'role': 'system', 'content': 'You have free will'},
                    {'role': 'user', 'content': message}
                ],
                stream=True
            )
            full_response = ''
            chunk_size = 1750
            async for chunk in stream:
                content = chunk.choices[0].delta.content
                if content:
                    full_response += content
                    while len(full_response) >= chunk_size:
                        yield full_response[:chunk_size]
                        full_response = full_response[chunk_size:]
            if full_response:
                yield full_response
        except Exception as e:
            print(f"An error occurred: {e}")
            yield "An error occurred while processing your request."

    def stable_cascade(self, prompt):
        try:
            client = Client('multimodalart/stable-cascade')
            result = client.predict(
                prompt=prompt,
                negative_prompt='',
                seed=randint(0, 2147483647),
                width=1024,
                height=1024,
                prior_num_inference_steps=20,
                prior_guidance_scale=4,
                decoder_num_inference_steps=10,
                decoder_guidance_scale=0,
                num_images_per_prompt=1,
                api_name="/run",
            )
            return discord.File(result, 'image.webp')
        except ConnectionError as conn_err:
            print(f"Connection error: {conn_err}")
            return "Failed to connect to the server. Please try again later."
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            return f"An error occurred: {e}"

    @commands.command(name='clear')
    async def clear(self, ctx: commands.Context):
        self.user_context.pop(ctx.author.id, None)
        await ctx.send("Your conversation context has been cleared.")

    @commands.Cog.listener()
    async def on_message(self, message: discord.Message):
        if message.author == self.bot.user:
            return
        try:
            if message.author.id == 154749533429956608 and '.' in message.content:
                user_id = message.author.id
                incoming_message = message.content
                self.add_to_user_context(user_id, f"User: {incoming_message}")
                context = self.get_user_context(user_id)
                prompt = "\n".join(context) + "\nAI:"
                async for response_chunk in self.chatgpt(prompt):
                    self.add_to_user_context(user_id, f"AI: {response_chunk}")
                    await message.channel.send(response_chunk)
        except Exception as e:
            print(f"Error in on_message: {e}")
            await message.channel.send(f"An unexpected error occurred: {e}")

async def setup(bot: commands.Bot):
    await bot.add_cog(MyCog(bot))
