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
from bot.main import CustomBot
from discord.ext import commands
from gradio_client import Client
from openai import OpenAI
from random import randint
import os

import discord

class MyCog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot

    async def openai(self, ctx: commands.Context):
        openai = OpenAI(api_key=self.bot.load_config()['api_keys']['api_key_1'])
        stream = self.openai.chat.completions.create(
             model='gpt-4o-mini',
             messages=[{'role': 'system', 'content': 'You are a chemistry expert.'},{'role': 'user', 'content': message.content}],
             stream=True,
             max_tokens=200
        )
        full_response = ''
        for chunk in stream:
            content = chunk.choices[0].delta.content
            if content:
                full_response += content
        return full_response

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

    @commands.Cog.listener()
    async def on_message(self, message: discord.Message):
        ctx = await self.bot.get_context(message)
        if message.author.bot:
            return
        #async with ctx.typing():
            #await ctx.send(await self.openai(ctx))
        try:
            file = self.stable_cascade(message.content)
            if isinstance(file, discord.File):
                await ctx.send(file=file)
            else:
                await ctx.send(f"Error generating image: {file}")
        except Exception as e:
            print(f"Error in on_message: {e}")
            await ctx.send(f"An unexpected error occurred: {e}")

async def setup(bot: commands.Bot):
    await bot.add_cog(MyCog(bot))
