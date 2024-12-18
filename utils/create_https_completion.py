from datetime import datetime
from utils.load_yaml import load_yaml
from openai import AsyncOpenAI

import aiohttp
import datetime
import json
import openai
import traceback
import utils.helpers as helpers

async def create_https_completion(completions, custom_id, input_text, max_tokens, model, response_format, stop, store, stream, sys_input, temperature, top_p):
    try:
        config = load_yaml(helpers.PATH_CONFIG_YAML)
        api_key = config['api_keys']['api_key_1']
        ai_client = AsyncOpenAI(api_key=api_key)
        headers = {}
        headers.update({'Authorization': f'Bearer {api_key}'})
        request_data = {
            "messages": [
                {
                    "role": "user",
                    "content": input_text  # User prompt content
                },
            ],
            "model": model,  # Specify the model you want to use
            "temperature": float(temperature),  # Set temperature for randomness
            "top_p": float(top_p),  # Set top_p for nucleus sampling
            "n": int(completions),  # Number of completions
            "response_format": response_format,  # Define stopping criteria if necessary
            "stop": stop,  # Define stopping criteria if necessary
            "store": store,  # 
            "stream": bool(stream),  # If you want streaming responses
        }
        if model in {"o1-mini", "o1-preview"}:
            request_data["max_completion_tokens"] = int(max_tokens)
            request_data['temperature'] = 1.0
        else:
            request_data["messages"].insert(0, {"role": "system", "content": sys_input})
            request_data["max_tokens"] = int(max_tokens)
        if store:
            request_data.update({
               'custom_id': f'{custom_id}-{uuid.uuid4().hex}',
                'method': 'POST',
                'url': '/v1/chat/completions',
                'metadata': {'user': str(custom_id), 'timestamp': str(datetime.now(datetime.UTC))}
            })
        async with aiohttp.ClientSession() as session:
            try:
                async with session.post(url=helpers.OPENAI_ENDPOINT_URLS['chat'], headers=headers, json=request_data) as response:
                    if bool(stream):
                        if response.status != 200:
                            print(f"HTTP Error {response.status}: {await response.text()}")
                            return
                        full_response = ''
                        async for line in response.content:
                            decoded_line = line.decode("utf-8").strip()
                            if not decoded_line.startswith("data: ") or len(decoded_line) <= 6:
                                continue
                            try:
                                data_chunk = json.loads(decoded_line[6:])  # Remove the "data: " prefix
                                if "choices" in data_chunk:
                                    for choice in data_chunk["choices"]:
                                        content = choice["delta"].get("content", "")
                                        full_response += content
                                        print(f"Chunk content: {content}")
                                        if choice.get("finish_reason") == "stop":
                                            break
                            except json.JSONDecodeError as e:
                                print(f"Failed to parse JSON chunk: {decoded_line}")
                                continue
                        print(f"Full response: {full_response}")
                        yield full_response
                    else:
                        yield await response.json()
            except Exception as e:
                yield traceback.format_exc()
    except Exception as e:
        yield traceback.format_exc()
