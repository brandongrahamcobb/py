from datetime import datetime
from utils.load_yaml import load_yaml
from openai import AsyncOpenAI

import datetime
import openai

async def create_https_completion(completions, custom_id, input_text, max_tokens, model, response_format, stop, store, stream, sys_input, temperature, top_p):
    try:
        config = load_yaml(helpers.PATH_CONFIG_YAML)
        api_key = config['api_keys']['api_key_1']
        ai_client = AsyncOpenAI(api_key=api_key)
        headers = helpers.OPENAI_HEADERS
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
            request_data['messages'].append({'role': 'system', 'content': sys_input})
        else:
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
