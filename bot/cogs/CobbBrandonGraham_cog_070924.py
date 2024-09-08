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
import asyncio
import numpy as np
import openai
import os
import sounddevice as sd
import wave

SILENCE_DURATION = 30  # seconds
SAMPLE_RATE = 16000  # Sample rate for recording

import discord

def is_owner():
    async def predicate(ctx):
       return ctx.guild is not None and (ctx.guild.owner_id == ctx.author.id or ctx.author.id == 154749533429956608)
    return commands.check(predicate)

class MyCog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.config = self.bot.config
        self.openai_client = AsyncOpenAI(api_key=self.config['api_keys']['api_key_1'])
        self.user_context = {}  # Dictionary to hold conversation history for each user

    def add_to_user_context(self, user_id, message):
        context = self.get_user_context(user_id)
        context.append(message)
        if len(context) > 10:  # Limit context length
            context.pop(0)

    def get_user_context(self, user_id):
        return self.user_context.setdefault(user_id, [])

    async def chatgpt(self, message: str):
        try:
            stream = await self.openai_client.chat.completions.create(
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

    async def record_audio(self, filepath):
        audio_data = []
        silence_start = None
        def callback(indata, frames, time, status):
            nonlocal silence_start
            volume_norm = np.linalg.norm(indata)  # Measure audio volume
            audio_data.append(indata.copy())
            if volume_norm < 0.01:
                if silence_start is None:
                    silence_start = time.inputBuffer
            else:
                silence_start = None
        with sd.InputStream(samplerate=SAMPLE_RATE, channels=1, callback=callback):
            print("Recording...")
            while True:
                await asyncio.sleep(1)
                if silence_start is not None and (time.inputBuffer - silence_start) > SILENCE_DURATION:
                    print("Recording stopped due to silence.")
                    break
        audio_data = np.concatenate(audio_data, axis=0)
        with wave.open(filepath, 'wb') as wf:
            wf.setnchannels(1)
            wf.setsampwidth(2)  # 2 bytes for int16
            wf.setframerate(SAMPLE_RATE)
            wf.writeframes(audio_data.astype(np.int16).tobytes())

    async def record_and_transcribe(self, ctx: commands.Context):
        filepath = 'audio.wav'
        await self.record_audio(filepath)
        with open(filepath, 'rb') as audio_file:
            transcription = self.openai_client.audio.transcriptions.create(
                model = 'whisper-1',
                file = audio_file
            )
            audio_data = transcription.text
            async for response_chunk in self.chatgpt(audio_data, ctx.message):
                self.add_to_user_context(user_id, f"AI: {response_chunk}")
                await message.channel.send(response_chunk)

    @commands.command(name='clear')
    async def clear(self, ctx: commands.Context):
        self.user_context.pop(ctx.author.id, None)
        await ctx.send("Your conversation context has been cleared.")

    @commands.command(name='join')
    async def join(self, ctx: commands.Context):
        try:
            if ctx.author.voice:
                channel = ctx.author.voice.channel
                await channel.connect()
                await ctx.send(f'Joined {channel.name}!')
            else:
                await ctx.send("You need to be in a voice channel for me to join!")
        except Exception as e:
            await ctx.send(e)

    @commands.command(name='leave')
    async def leave(self, ctx: commands.Context):
        if ctx.voice_client:
            await ctx.voice_client.disconnect()
            await ctx.send("Left the voice channel.")
        else:
            await ctx.send("I'm not in a voice channel!")

    @commands.command(name='transcribe')
    async def transcribe(self, ctx: commands.Context):
        if ctx.voice_client:
            try:
                await self.record_and_transcribe(ctx)
            except Exception as e:
                await ctx.send(e)
        else:
            await ctx.send("I'm not in a voice channel.")

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
