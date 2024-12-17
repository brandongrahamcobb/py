from openai import AsyncOpenAI

import load_yaml
import openai
import traceback

async def create_completion(input_text, conversation_id, conversations, model):
    try:
        config = load_yaml(helpers.PATH_CONFIG_YAML)
        api_key = config['api_keys']['api_key_1']
        ai_client = AsyncOpenAI(api_key=api_key)
        messages = conversations[conversation_id]
        messages.append({'role': 'users', 'content': input_text})
        stream = await ai_client.chat.completions.create(
            model=model,
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
