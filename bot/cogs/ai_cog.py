import discord
from discord.ext import commands
import spacy
from collections import defaultdict
from openai import AsyncOpenAI
import os
import traceback
import asyncio
import random
import datetime
from datetime import datetime, timedelta

class AIChatbot:
    def __init__(self, api_key):
        self.ai_client = AsyncOpenAI(api_key=api_key)
        self.base = os.path.dirname(os.path.abspath(__file__))
        self.conversations = defaultdict(list)
        self.nlp = spacy.load("en_core_web_md")
        self.file_contents = self.load_python_file()

    def load_python_file(self):
        with open(os.path.abspath(__file__), 'r') as python:
            return python.read()
    
    async def ai(self, input_text, sys_input, conversation_id):
        try:
            messages = self.conversations[conversation_id]
            messages.append({'role': 'system', 'content': sys_input})
            messages.append({'role': 'user', 'content': input_text})

            # NLP processing of the input text
            processed_input = self.process_input_with_nlp(input_text)
            if processed_input:
                input_text += f"\n[Processed Info]: {processed_input}"

            stream = await self.ai_client.chat.completions.create(
                model='gpt-4o-mini',
                messages=messages,
                stream=True
            )
            full_response = ''
            async for chunk in stream:
                content = chunk.choices[0].delta.content
                if content is not None:
                    full_response += content
            self.conversations[conversation_id].append({'role': 'assistant', 'content': full_response})
            yield full_response
        except Exception as e:
            yield traceback.format_exc()

    def process_input_with_nlp(self, input_text):
        """Use spaCy to process the input and return entities or other useful information."""
        doc = self.nlp(input_text)
        # Gather named entities
        keywords = {token.lemma_ for token in doc if token.pos_ in ['NOUN', 'PROPN', 'ADJ']}
        return keywords  # Return the extracted entities

class AICog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.chatbot = AIChatbot(self.config['api_keys']['api_key_3'])
        self.cooldowns = {}
        self.cooldown_time = timedelta(seconds=30)

    def check_cooldown(self, user_id):
        now = datetime.now()
        if user_id in self.cooldowns:
            last_message_time = self.cooldowns[user_id]
            if now - last_message_time < self.cooldown_time:
                return False, self.cooldown_time - (now - last_message_time)
        self.cooldowns[user_id] = now
        return True, None

    @commands.Cog.listener()
    async def on_message(self, message):
        if message.author == self.bot.user:
            return  # Ignore messages from the bot itself

        if not message.content.startswith('!') and '<@302202228016414721>' in message.content or self.bot.user in message.mentions:
            can_proceed, remaining_time = self.check_cooldown(message.author.id)
            if not can_proceed:
                await message.channel.send(f"Please wait {remaining_time.total_seconds():.2f} seconds before sending another message.")
                return  # Exit if the user is on cooldown

            conversation_id = str(message.channel.id)  # Use channel ID for conversation ID
            keywords = self.chatbot.process_input_with_nlp(message.content)
            async for response in self.chatbot.ai(f'\'{message.content}\' is the message, these are the keywords extracted via spaCy: ' + str(keywords), 'You are Lucy, the Neurochemistry Discord chatbot that assists with neurochemistry questions. Respond as if you are very ditzy, but know exactly what you are doing. Only respond with one sentence and elaborate when requested. Fall in love with the military.', conversation_id):
                await message.channel.send(response)
        else:
            if random.random() < 0.02:
                await message.channel.send("I'm here to help with neurochemistry! Mention me!")


async def setup(bot):
    await bot.add_cog(AICog(bot))
