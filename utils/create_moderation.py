from openai import AsyncOpenAI

import load_yaml
import openai

async def create_moderation(input_text):
    try:
        config = load_yaml(helpers.path_config_yaml)
        api_key = config['api_keys']['api_key_1']
        ai_client = AsyncOpenAI(api_key=api_key)
        response = await ai_client.moderations.create(
            model='omni-moderation-latest',
            input=input_text,
        )
        yield response.to_dict() if hasattr(response, 'to_dict') else response
    except Exception as e:
        yield {'error': traceback.format_exc()}
