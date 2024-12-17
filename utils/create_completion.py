from openai import AsyncOpenAI

import load_yaml
import openai

async def create_completion_deprecated(input_text, conversation_id, model, sys_input):
    try:
        config = helpers.load_yaml(helpers.path_config_yaml)
        api_key = config['api_keys']['api_key_1']
        ai_client = AsyncOpenAI(api_key=api_key)
        messages = conversations[conversation_id]
        messages.append({'role': 'system', 'content': sys_input})
        messages.append({'role': 'user', 'content': input_text})
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
