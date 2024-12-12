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

conversations = defaultdict(list)
home = expanduser('~')
path_config_yaml = join(home, '.config', 'vyrtuous', 'config.yaml')

def api_endpoint_from_url(request_url):
    match = re.search("^https://[^/]+/v\\d+/(.+)$", request_url)
    return match[1]

def append_to_jsonl(data, filename: str) -> None:
    json_string = json.dumps(data)
    with open(filename, "a") as f:
        f.write(json_string + "\n")

@dataclass
class StatusTracker:
    num_tasks_started: int = 0
    num_tasks_in_progress: int = 0
    num_tasks_succeeded: int = 0
    num_tasks_failed: int = 0
    num_rate_limit_errors: int = 0
    num_api_errors: int = 0
    num_other_errors: int = 0
    time_of_last_rate_limit_error: float = 0.0


@dataclass
class APIRequest:
    task_id: int
    request_json: dict
    token_consumption: int
    attempts_left: int
    metadata: dict
    result: list = field(default_factory=list)
    async def call_api(
        self,
        session: aiohttp.ClientSession,
        request_url: str,
        request_header: dict,
        retry_queue: asyncio.Queue,
        save_filepath: str,
        status_tracker: StatusTracker,
    ):
        logging.info(f"Starting request #{self.task_id}")
        error = None
        try:
            async with session.post(
                url=request_url, headers=request_header, json=self.request_json
            ) as response:
                response = await response.json()
            if "error" in response:
                logging.warning(f"Request {self.task_id} failed with error {response['error']}")
                status_tracker.num_api_errors += 1
                error = response
                if "rate limit" in response["error"].get("message", "").lower():
                    status_tracker.time_of_last_rate_limit_error = time.time()
                    status_tracker.num_rate_limit_errors += 1
                    status_tracker.num_api_errors -= 1  
        except Exception as e:
            logging.warning(f"Request {self.task_id} failed with Exception {e}")
            status_tracker.num_other_errors += 1
            error = e
        if error:
            self.result.append(error)
            if self.attempts_left:
                retry_queue.put_nowait(self)
            else:
                logging.error(f"Request {self.request_json} failed after all attempts.")
                data = [self.request_json, [str(e) for e in self.result], self.metadata] if self.metadata else [self.request_json, [str(e) for e in self.result]]
                append_to_jsonl(data, save_filepath)
                status_tracker.num_tasks_in_progress -= 1
                status_tracker.num_tasks_failed += 1
        else:
            data = [self.request_json, response, self.metadata] if self.metadata else [self.request_json, response]
            append_to_jsonl(data, save_filepath)
            status_tracker.num_tasks_in_progress -= 1
            status_tracker.num_tasks_succeeded += 1
            logging.debug(f"Request {self.task_id} saved to {save_filepath}")

async def create_completion(completions, conversation_id, input_text, max_tokens, model, stop, store, stream, sys_input, temperature, top_p):
    try:
        config = helpers.load_yaml(path_config_yaml)
        api_key = config['api_keys']['api_key_2']
        ai_client = AsyncOpenAI(api_key=api_key)
        url = "https://api.openai.com/v1/chat/completions"
        headers = {
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json",
            "OpenAI-Organization": "spawd",
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
        if store:
            request_data.update({
                "metadata": {
                    "user": str(conversation_id),
                    "timestamp": str(datetime.utcnow())
                }
            })
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
        config = helpers.load_yaml(path_config_yaml)
        api_key = config['api_keys']['api_key_2']
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
        config = helpers.load_yaml(path_config_yaml)
        api_key = config['api_keys']['api_key_2']
        ai_client = AsyncOpenAI(api_key=api_key)
        messages = conversations[conversation_id]
        messages.append({'role': 'system', 'content': sys_input})
        messages.append({'role': 'user', 'content': input_text})
        stream = await ai_client.chat.completions.create(
            model='gpt-3.5-turbo',
            messages=messages,
            stream=True
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
        config = helpers.load_yaml(path_config_yaml)
        api_key = config['api_keys']['api_key_2']
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

def num_tokens_consumed_from_request(request_json: dict, api_endpoint: str, token_encoding_name: str):
    encoding = tiktoken.get_encoding(token_encoding_name)
    if api_endpoint.endswith("completions"):
        max_tokens = request_json.get("max_tokens", 15)
        n = request_json.get("n", 1)
        completion_tokens = n * max_tokens
        prompt = request_json["prompt"]
        if isinstance(prompt, str):
            prompt_tokens = len(encoding.encode(prompt))
            return prompt_tokens + completion_tokens

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

async def process_api_requests_from_file(
    requests_filepath: str,
    save_filepath: str,
    request_url: str,
    max_requests_per_minute: float,
    max_tokens_per_minute: float,
    token_encoding_name: str,
    max_attempts: int,
    logging_level: int,
):
    """Processes API requests in parallel, throttling to stay under rate limits."""
    seconds_to_pause_after_rate_limit_error = 15
    seconds_to_sleep_each_loop = 0.001  
    logging.basicConfig(level=logging_level)
    logging.debug(f"Logging initialized at level {logging_level}")
    api_endpoint = api_endpoint_from_url(request_url)
    request_header = {"Authorization": f"Bearer {config['api_keys']['api_key_2']}"}
    queue_of_requests_to_retry = asyncio.Queue()
    task_id_generator = task_id_generator_function()
    status_tracker = StatusTracker()
    next_request = None  
    available_request_capacity = max_requests_per_minute
    available_token_capacity = max_tokens_per_minute
    last_update_time = time.time()
    file_not_finished = True  
    with open(requests_filepath) as file:
        requests = file.__iter__()
        logging.debug(f"File opened. Entering main loop")
        async with aiohttp.ClientSession() as session:
            while True:
                if next_request is None:
                    if not queue_of_requests_to_retry.empty():
                        next_request = queue_of_requests_to_retry.get_nowait()
                    elif file_not_finished:
                        try:
                            request_json = json.loads(next(requests))
                            next_request = APIRequest(
                                task_id=next(task_id_generator),
                                request_json=request_json,
                                token_consumption=num_tokens_consumed_from_request(
                                    request_json, api_endpoint, token_encoding_name
                                ),
                                attempts_left=max_attempts,
                                metadata=request_json.pop("metadata", None),
                            )
                            status_tracker.num_tasks_started += 1
                            status_tracker.num_tasks_in_progress += 1
                        except StopIteration:
                            file_not_finished = False
                current_time = time.time()
                seconds_since_update = current_time - last_update_time
                available_request_capacity = min(
                    available_request_capacity + max_requests_per_minute * seconds_since_update / 60.0,
                    max_requests_per_minute,
                )
                available_token_capacity = min(
                    available_token_capacity + max_tokens_per_minute * seconds_since_update / 60.0,
                    max_tokens_per_minute,
                )
                last_update_time = current_time
                if next_request:
                    next_request_tokens = next_request.token_consumption
                    if available_request_capacity >= 1 and available_token_capacity >= next_request_tokens:
                        available_request_capacity -= 1
                        available_token_capacity -= next_request_tokens
                        next_request.attempts_left -= 1
                        asyncio.create_task(
                            next_request.call_api(
                                session=session,
                                request_url=request_url,
                                request_header=request_header,
                                retry_queue=queue_of_requests_to_retry,
                                save_filepath=save_filepath,
                                status_tracker=status_tracker,
                            )
                        )
                        next_request = None  
                if status_tracker.num_tasks_in_progress == 0:
                    break

                await asyncio.sleep(seconds_to_sleep_each_loop)
                seconds_since_rate_limit_error = (
                    time.time() - status_tracker.time_of_last_rate_limit_error
                )
                if (
                    seconds_since_rate_limit_error
                    < seconds_to_pause_after_rate_limit_error
                ):
                    remaining_seconds_to_pause = (
                        seconds_to_pause_after_rate_limit_error
                        - seconds_since_rate_limit_error
                    )
                    await asyncio.sleep(remaining_seconds_to_pause)
                    logging.warn(
                        f"Pausing to cool down until {time.ctime(status_tracker.time_of_last_rate_limit_error + seconds_to_pause_after_rate_limit_error)}"
                    )
        logging.info(f"Parallel processing complete. Results saved to {save_filepath}")
        if status_tracker.num_tasks_failed > 0:
            logging.warning(
                f"{status_tracker.num_tasks_failed} / {status_tracker.num_tasks_started} requests failed."
            )
        if status_tracker.num_rate_limit_errors > 0:
            logging.warning(
                f"{status_tracker.num_rate_limit_errors} rate limit errors received."
            )


