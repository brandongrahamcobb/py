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

import aiohttp
import datetime
import json
import openai
import traceback
import utils.helpers as helpers

async def create_https_completion(completions, conversations, custom_id, input_array, max_tokens, model, response_format, stop, store, stream, sys_input, temperature, top_p):
    try:
        config = load_yaml(helpers.PATH_CONFIG_YAML)
        api_key = config['api_keys']['api_key_1']
        ai_client = AsyncOpenAI(api_key=api_key)
        headers = {}
        headers.update({'Authorization': f'Bearer {api_key}'})
        request_data = {
            'messages': conversations[custom_id],
            'model': model,
            'temperature': float(temperature),
            'top_p': float(top_p),
            'n': int(completions),
            'stop': stop,
            'store': bool(store),
            'stream': bool(stream),
        }
        last = len(request_data['messages']) - 1
        request_data['messages'].insert(last, {'role': 'user', 'content': input_array})
        if response_format:
            request_data['response_format'] = response_format
        if model in {'o1-mini', 'o1-preview'}:
            request_data['max_completion_tokens'] = int(max_tokens)
            request_data['temperature'] = 1.0
        else:
            request_data['messages'].insert(0, {'role': 'system', 'content': sys_input})
            request_data['max_tokens'] = int(max_tokens)
        if bool(store):
            request_data.update({
#               'custom_id': f'{custom_id}-{uuid.uuid4().hex}',
 #               'method': 'POST',
  #              'url': '/v1/chat/completions',
                'metadata': {'user': str(custom_id), 'timestamp': str(datetime.datetime.now(datetime.timezone.utc))}
            })
        async with aiohttp.ClientSession() as session:
            try:
                async with session.post(url=helpers.OPENAI_ENDPOINT_URLS['chat'], headers=headers, json=request_data) as response:
                    if bool(stream):
                        if response.status != 200:
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
                                            break
                            except json.JSONDecodeError as e:
                                continue
                        yield full_response
                    else:
                        yield await response.json()
            except Exception as e:
                yield traceback.format_exc()
    except Exception as e:
        yield traceback.format_exc()

