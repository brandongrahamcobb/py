''' openai_helpers.py
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
'''
from collections import defaultdict
from datetime import datetime
from dataclasses import (
    dataclass,
    field,
)
from openai import AsyncOpenAI
from typing import Any, List, Optional, Dict

import aiohttp
import asyncio
import bot.utils.helpers as helpers
import discord
import json
import logging
import re
import tiktoken
import time
import traceback
from os.path import expanduser, join
import os
import uuid

conversations = defaultdict(list)
home = expanduser('~')

async def create_completion(completions, conversation_id, input_text, max_tokens, model, stop, store, stream, sys_input, temperature, top_p):
    try:
        config = helpers.load_yaml(helpers.path_config_yaml)
        api_key = config['api_keys']['api_key_1']
        ai_client = AsyncOpenAI(api_key=api_key)
        url = "https://api.openai.com/v1/chat/completions"
        headers = {
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json",
            "OpenAI-Organization": "org-3LYwtg7DSFJ7RLn9bfk4hATf",
            "OpenAI-Project": "proj_u5htBCWX0LSHxkw45po1Vfz9",
        }
        request_data = {
            "max_tokens": int(max_tokens),  # Set maximum tokens
            "messages": [
                {
                    "role": "user",
                    "content": input_text  # User prompt content
                },
                {
                    "role": "system",
                    "content": sys_input  # User prompt content
                },
            ],
            "model": model,  # Specify the model you want to use
            "temperature": float(temperature),  # Set temperature for randomness
            "top_p": float(top_p),  # Set top_p for nucleus sampling
            "n": int(completions),  # Number of completions
            "stop": stop,  # Define stopping criteria if necessary
            "store": store,  # 
            "stream": bool(stream),  # If you want streaming responses
        }
        if model in {"o1-mini", "o1-preview"}:
            request_data["max_completion_tokens"] = int(max_tokens)
        else:
            request_data["max_tokens"] = int(max_tokens)
        if store:
            return
        async with aiohttp.ClientSession() as session:
            try:
                async with session.post(url, headers=headers, json=request_data) as response:
                    if response.status != 200:
                        logging.error(f"HTTP Error: {response.status} - {await response.text()}")
                        yield f"HTTP error {response.status}: {await response.text()}"
                        return
                    if stream:
                        async for line in response.content:
                            decoded_line = line.decode("utf-8").strip()
                            if decoded_line.startswith("data: "):
                                data = decoded_line[6:]
                                if data == "[DONE]":
                                    break
                                try:
                                    json_data = json.loads(data)
                                    full_response = ''
                                    if "choices" in json_data:
                                        for choice in json_data["choices"]:
                                            content = choice["delta"].get("content")
                                            if content:
                                                full_response += content
                                            else:
                                                logging.warning("Missing 'content' in the response.")
                                            yield full_response
                                except json.JSONDecodeError as e:
                                    logging.error(f"JSON decode error: {e} - Line: {decoded_line}")
                    else:
                        result = await response.json()
                        if "choices" in result:
                            for choice in result["choices"]:
                                yield choice.get("message", {}).get("content", "")
            except Exception as e:
                yield traceback.format_exc()
    except Exception as e:
        yield traceback.format_exc()

async def create_moderation(input_text):
    try:
        config = helpers.load_yaml(helpers.path_config_yaml)
        api_key = config['api_keys']['api_key_1']
        ai_client = AsyncOpenAI(api_key=api_key)
        response = await ai_client.moderations.create(
            model='omni-moderation-latest',
            input=input_text,
        )
        yield response.to_dict() if hasattr(response, 'to_dict') else response
    except Exception as e:
        yield {'error': traceback.format_exc()}

async def create_pseudomoderation(input_text, sys_input, conversation_id):
    try:
        config = helpers.load_yaml(helpers.path_config_yaml)
        api_key = config['api_keys']['api_key_1']
        ai_client = AsyncOpenAI(api_key=api_key)
        messages = conversations[conversation_id]
        messages.append({'role': 'system', 'content': sys_input})
        messages.append({'role': 'user', 'content': input_text})
        stream = await ai_client.chat.completions.create(
            model='gpt-3.5-turbo',
            messages=messages,
            stream=True,
            response_format={ "type": "json_object" }
        )
        full_response = ''
        async for chunk in stream:
            content = chunk.choices[0].delta.content
            if content is not None:
                full_response += content
        conversations[conversation_id].append({'role': 'assistant', 'content': full_response})
        yield full_response
    except Exception as e:
        yield traceback.format_exc()

async def create_transcription(audio_bytes, sample_rate=48000):
    try:
        config = helpers.load_yaml(helpers.path_config_yaml)
        api_key = config['api_keys']['api_key_1']
        ai_client = AsyncOpenAI(api_key=api_key)
        audio = AudioSegment(
            data=audio_bytes,
            sample_width=2,          # 16-bit audio
            frame_rate=sample_rate,  # Typically 48000 for Discord
            channels=1                # Mono
        )
        with NamedTemporaryFile(suffix=".wav") as temp_audio_file:
            audio.export(temp_audio_file.name, format="wav")
            with open(temp_audio_file.name, "rb") as f:
                transcript = ai_client.Audio.transcribe(
                    model="whisper-1",
                    file=f
                )
        return transcript.get("text", "")
    except Exception as e:
        print(f"Error during transcription: {e}")
        return "Transcription failed."

async def handle_module_moderation(message: discord.Message, moderation):
    if isinstance(moderation, dict) and 'error' in moderation:
        print(f"Moderation error: {moderation['error']}")
        return  # Handle errors (e.g., log or notify users)
     # Check the structure of the moderation response
    if 'results' in moderation and len(moderation['results']) > 0:
        flagged = moderation['results'][0].get('flagged', False)
        print("Flagged Status:", flagged)
        categories = moderation['results'][0].get('categories', {})
        print("Categories:", categories)
        if bool(flagged) or flagged:
            await message.delete()  # Delete the flagged message
            await message.channel.send(f"{message.author.mention}'s post was flagged and deleted due to inappropriate content.")

async def handle_api_moderation(message: discord.Message, moderation):
    if not moderation or not moderation.strip():
        print("Received an empty or null response from moderation.")
    moderation_dict = {}
    try:
        moderation_dict = json.loads(moderation)
    except json.JSONDecodeError as e:
        print(f"JSON decoding failed: {e}")
        print(f"Response content was: '{moderation}")
    flagged = moderation_dict.get('results', [{}])[0].get('flagged', False)
    if bool(flagged) or flagged:
        spoiler_content = f"||{message.content}||"
        await message.delete()
        await message.channel.send(f"{message.author.mention} said: {spoiler_content}")


async def process_batch_from_file(batch_input_filename="requests.jsonl"):
    """
    Processes a batch using OpenAI's Batch API from the given JSONL file.
    
    Args:
        batch_input_filename (str): Path to the requests.jsonl file.
    
    Returns:
        tuple: (success (bool), message (str)) where success is True on success, and message contains the results or error details.
    """
    config = helpers.load_yaml(helpers.path_config_yaml)
    api_key = config['api_keys']['api_key_1']
    ai_client = AsyncOpenAI(api_key=api_key)
    try:
        # Ensure the file exists and load requests
        with open(batch_input_filename, "r") as file:
            requests = [json.loads(line) for line in file]
    except FileNotFoundError:
        return False, f"File {batch_input_filename} not found. Please create the file and try again."
    except json.JSONDecodeError:
        return False, f"Error decoding {batch_input_filename}. Ensure it's in proper JSONL format."

    # Step 1: Upload the batch file
    try:
        batch_input_file = await ai_client.files.create(file=open(batch_input_filename, "rb"), purpose="batch")
        batch_input_file_id = batch_input_file.id
    except Exception as e:
        return False, f"Error uploading batch file: {e}"

    # Step 2: Create the batch
    try:
        batch = await ai_client.batches.create(
            input_file_id=batch_input_file_id,
            endpoint="/v1/chat/completions",
            completion_window="24h",
            metadata={"description": "Discord batch processing"}
        )
        batch_id = batch.id  # Use attribute-style access here
    except Exception as e:
        return False, f"Error creating batch: {e}"

    # Step 3: Monitor the batch status
    while True:
        batch_status = await ai_client.batches.retrieve(batch_id)
        print(batch_status.errors)
        if batch_status.status == "completed":
            break
        elif batch_status.status in ["failed", "cancelled", "expired"]:
            return False, f"Batch failed with status: {batch_status.status}"

    # Step 4: Retrieve the results
    try:
        output_file_id = batch_status.output_file_id  # Use attribute-style access
        output_file = await ai_client.files.content(output_file_id)
        output_data = [json.loads(line) for line in output_file.splitlines()]

        # Generate results string
        results = "\n".join(
            f"Request {i + 1}: {res['response']['body']['choices'][0]['message']['content']}"
            for i, res in enumerate(output_data)
        )
        return True, results
    except Exception as e:
        return False, f"Error retrieving batch results: {e}"


def get_user_preferences(user_id: str):
    """Fetch user-specific preferences from storage."""
    prefs = helpers.load_json(helpers.path_preferences) or {}
    return prefs.get(user_id, {
        "completions": 1,
        "conversation_id": "",
        "max_tokens": 150,
        "model": "gpt-3.5-turbo",
        "stop": None,
        "store": False,
        "stream": True,
        "sys_input": "You are a helpful assistant.",
        "temperature": 0.7,
        "top_p": 1.0,
    })

# Embed Templates for Tutorial Sections
def get_prerequisites_embed():
    embed = discord.Embed(
        title="Step 1: Prerequisites",
        description="Before we start, make sure you have the following:",
        color=discord.Color.blue(),
    )
    embed.add_field(name="1. Python", value="Download and install Python 3.8+ from [python.org](https://www.python.org/).", inline=False)
    embed.add_field(name="2. Install Required Libraries", value="Run this command: `pip install discord.py openai`", inline=False)
    embed.set_footer(text="Use the dropdown below to navigate the tutorial!")
    return embed

def get_api_key_embed():
    embed = discord.Embed(
        title="Step 2: Get Your OpenAI API Key",
        description=(
            "1. Go to [OpenAI](https://platform.openai.com/).\n"
            "2. Sign up or log in.\n"
            "3. Navigate to the API keys section in your account.\n"
            "4. Copy your API key and save it securely."
        ),
        color=discord.Color.green(),
    )
    embed.set_footer(text="Use the dropdown below to navigate the tutorial!")
    return embed

def get_preferences_embed(user_id: str):
    preferences = get_user_preferences(user_id)
    embed = discord.Embed(
        title="Your Current Settings",
        description="These are your current settings for OpenAI completions:",
        color=discord.Color.purple(),
    )
    for key, value in preferences.items():
        embed.add_field(name=key.capitalize(), value=str(value), inline=False)
    embed.set_footer(text="Use the dropdown below to navigate the tutorial!")
    return embed

# Tutorial View with Dropdown
class TutorialDropdown(discord.ui.Select):
    def __init__(self, user_id: str):
        self.user_id = user_id
        options = [
            discord.SelectOption(label="Step 1: Prerequisites", description="Learn what you need before starting.", emoji="📋"),
            discord.SelectOption(label="Step 2: Get API Key", description="Learn how to get your OpenAI API key.", emoji="🔑"),
            discord.SelectOption(label="Your Current Settings", description="View your saved preferences.", emoji="⚙️"),
        ]
        super().__init__(placeholder="Select a step to learn...", max_values=1, min_values=1, options=options)

    async def callback(self, interaction: discord.Interaction):
        if self.values[0] == "Step 1: Prerequisites":
            await interaction.response.edit_message(embed=get_prerequisites_embed())
        elif self.values[0] == "Step 2: Get API Key":
            await interaction.response.edit_message(embed=get_api_key_embed())
        elif self.values[0] == "Your Current Settings":
            await interaction.response.edit_message(embed=get_preferences_embed(self.user_id))

class TutorialView(discord.ui.View):
    def __init__(self, user_id: str):
        super().__init__()
        self.add_item(TutorialDropdown(user_id))
