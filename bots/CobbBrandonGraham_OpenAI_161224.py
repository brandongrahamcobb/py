from collections import defaultdict
from datetime import datetime as dt
from os.path import join

import aiohttp
import api_request_parallel_processor as arpp
import asyncio
import helpers
import json
import traceback
import uuid

CONFIG = helpers.load_yaml(helpers.path_config_yaml)
API_KEY = CONFIG['api_keys']['api_key_1']
conversations = defaultdict(list)
COMPLETIONS = 1
CUSTOM_ID = uuid.uuid4().hex
HEADERS = {
    'Authorization': f'Bearer {API_KEY}',
    'Content-Type': 'application/json',
    'OpenAI-Project': 'proj_e3d8jKn83rRCHhf2CZi1xv9a',
    'User-Agent': 'brandongrahamcobb@icloud.com'
}
INPUT_TEXT = input('INPUT_TEXT = ')  # Ensure that this line gets a string input

MODEL_OUTPUT_LIMITS = {
    "gpt-3.5-turbo": 4096,
    "gpt-4": 8192,
    "gpt-4-32k": 32768,
    "gpt-4o": 4096,         # Initially capped at 4,096; updated to 16,384 in later versions
    "gpt-4o-mini": 16384,
    "gpt-4-turbo": 4096,
    "o1-preview": 32768,
    "o1-mini": 65536,
}

MODEL_CONTEXT_LIMITS = {
    "gpt-3.5-turbo": 4096,
    "gpt-4": 8192,
    "gpt-4-32k": 32768,
    "gpt-4o": 128000,
    "gpt-4o-mini": 128000,
    "gpt-4-turbo": 128000,
    "o1-preview": 128000,
    "o1-mini": 128000,
}
MODELS = {
    'current': ['o1-preview', 'o1-mini'],
    'deprecated': ['gpt-3.5-turbo', 'gpt-4', 'gpt-4-32k', 'gpt-4o', 'gpt-4o-mini', 'gpt-4-turbo']
}
MODEL = 'o1-preview'
STOP = None
STORE = True
STREAM = False
SYS_INPUT = None  # Make sure this is defined correctly based on your usage
TOP_P = 1.0
TEMPERATURE = 0.7 if MODEL in MODELS['deprecated'] else 1.0

REQUEST_DATA = {
    "metadata": {
        "user": 'Brandon Graham Cobb',
        "timestamp": str(dt.utcnow())
    },
    "model": MODEL,
    "n": int(COMPLETIONS),
    "stop": STOP,
    "store": STORE,
    "stream": bool(STREAM),
    "temperature": float(TEMPERATURE),
    "top_p": float(TOP_P),
}

MAX_TOKENS = MODEL_OUTPUT_LIMITS[MODEL]

if MODEL in MODELS['current']:
    REQUEST_DATA.update({
        "max_completion_tokens": int(MAX_TOKENS),
        'messages': [{'role': 'user', 'content': INPUT_TEXT}]
    })
else:
    REQUEST_DATA.update({"max_tokens": int(MAX_TOKENS)})
    if SYS_INPUT:  # Include the system message only if it's defined
        REQUEST_DATA.update({
            'messages': [
                {'role': 'user', 'content': INPUT_TEXT},
                {'role': 'system', 'content': SYS_INPUT}
            ]
        })
    else:
        REQUEST_DATA.update({
            'messages': [{'role': 'user', 'content': INPUT_TEXT}]
        })


URL = "https://api.openai.com/v1/chat/completions"

async def main():
    async with aiohttp.ClientSession() as session:
        try:
            async with session.post(url=URL, headers=HEADERS, json=REQUEST_DATA) as response:
                if response.status == 200:
                    result = await response.json()
                    full_response = ''
                    if "choices" in result:
                        for choice in result["choices"]:
                            full_response = choice.get("message", {}).get("content", "")
                            conversations[CUSTOM_ID].append({'role': 'assistant', 'content': full_response})
                    return full_response
                else:
                    return f"Error: {response.status} - {await response.text()}"
        except Exception as e:
            return f"Exception: {str(e)}"

if __name__ == '__main__':
    response = asyncio.run(main())
    if response:  # Check if the response is not empty
        print(response.strip())
        path_responses = join(helpers.path_home, 'response.txt')
        with open(path_responses, 'a') as f:
            f.write(response.strip() + "\n")  # Append a newline for clarity
