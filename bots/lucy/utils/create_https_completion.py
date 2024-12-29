''' create_https_completion.py  OpenAI's v1/chat/completions endpoint using their python SDK is
                                much more efficient than this program. This is a complicated
                                way to get a completion from OpenAI from cd ../.
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

from collections import defaultdict
from datetime import datetime
from openai import AsyncOpenAI
from utils.load_yaml import load_yaml
from utils.setup_logging import logger

import aiohttp
import datetime
import json
import openai
import traceback
import utils.helpers as helpers

class Conversations:
    def __init__(self):
        self.conversations = defaultdict(list)

    def trim_conversation_history(self, model, custom_id):
        max_context_length = helpers.OPENAI_MODEL_CONTEXT_LIMITS.get(model, 4096)
        total_tokens = sum(len(msg['content']) for msg in self.conversations[custom_id])
        while total_tokens > max_context_length:
            removed_message = self.conversations[custom_id].pop(0)
            total_tokens -= len(removed_message['content'])
        
    async def create_https_completion(self, completions, custom_id, input_array, max_tokens, model, response_format, stop, store, stream, sys_input, temperature, top_p):
        try:
            logger.info("Loading configuration file.")
            config = load_yaml(helpers.PATH_CONFIG_YAML)
            api_key = config['api_keys']['api_key_1']['api_key']
            logger.info("API key loaded successfully.")
    
            ai_client = AsyncOpenAI(api_key=api_key)
            headers = {'Authorization': f'Bearer {api_key}'}
            logger.info("Headers prepared for the request.")
    
            request_data = {
                'messages': self.conversations[custom_id],
                'model': model,
                'temperature': float(temperature),
                'top_p': float(top_p),
                'n': int(completions),
                'stop': stop,
                'store': bool(store),
                'stream': bool(stream),
            }
            logger.info(f"Request data initialized for model: {model}.")
    
            last = len(request_data['messages']) - 1
            request_data['messages'].insert(last, {'role': 'user', 'content': input_array})
            logger.info("User input added to the messages.")
    
            if response_format:
                request_data['response_format'] = response_format
                logger.info(f"Response format set: {response_format}.")
    
            if model in {'chatgpt-4o-latest, o1-mini', 'o1-preview'}:
                request_data['max_completion_tokens'] = int(max_tokens)
                request_data['temperature'] = 1.0
                logger.info("Special settings applied for model: o1-mini or o1-preview.")
            else:
                request_data['messages'].insert(0, {'role': 'system', 'content': sys_input})
                request_data['max_tokens'] = int(max_tokens)
                logger.info("System input and max tokens added to the request data.")
    
            if bool(store):
                request_data.update({
                    'metadata': {'user': str(custom_id), 'timestamp': str(datetime.datetime.now(datetime.timezone.utc))}
                })
                logger.info("Store option enabled, metadata added to request data.")
    
            async with aiohttp.ClientSession() as session:
                try:
                    logger.info("Sending request to OpenAI chat endpoint.")
                    async with session.post(url=helpers.OPENAI_ENDPOINT_URLS['chat'], headers=headers, json=request_data) as response:
                        logger.info(f"Received response with status: {response.status}.")
    
                        if bool(stream):
                            if response.status != 200:
                                logger.error("Streaming response status not 200. Exiting.")
                                return
                            full_response = ''
                            async for line in response.content:
                                decoded_line = line.decode('utf-8').strip()
                                if not decoded_line.startswith('data: ') or len(decoded_line) <= 6:
                                    continue
                                try:
                                    data_chunk = json.loads(decoded_line[6:])  # Remove the 'data: ' prefix
                                    if 'choices' in data_chunk:
                                        for choice in data_chunk['choices']:
                                            content = choice['delta'].get('content', '')
                                            full_response += content
                                            if choice.get('finish_reason') == 'stop':
                                                logger.info("Completion streaming stopped.")
                                                break
                                except json.JSONDecodeError as e:
                                    logger.warning("Failed to decode JSON chunk during streaming.")
                                    continue
                            logger.info("Streaming response processed successfully.")
                        else:
                            logger.info("Processing non-streaming response.")
                            full_response_json = await response.json()
                            full_response = full_response_json['choices'][0]['message']['content']
                        if len(full_response) > helpers.DISCORD_CHARACTER_LIMIT:
                            char_limit = helpers.DISCORD_CHARACTER_LIMIT
                            is_code_block = False
                            parts = full_response.split("```")
                            for i in range(len(parts)):
                                if is_code_block:
                                    code_block_chunks = [parts[i][j:j+char_limit] for j in range(0, len(parts[i]), char_limit)]
                                    for chunk in code_block_chunks:
                                        if self.is_replying_all == "True" or has_followed_up:
                                            yield f'```{chunk}```'
                                        else:
                                            yield f'```{chunk}```'
                                    is_code_block = False
                                else:
                                    non_code_chunks = [parts[i][j:j+char_limit] for j in range(0, len(parts[i]), char_limit)]
                                    for chunk in non_code_chunks:
                                       yield chunk
                        else:
                            yield full_response
                except Exception as e:
                    logger.error("An error occurred while making the HTTP request.", exc_info=True)
                    yield traceback.format_exc()
        except Exception as e:
            logger.error("An error occurred in create_https_completion.", exc_info=True)
            yield traceback.format_exc()
        self.trim_conversation_history(model, custom_id)
