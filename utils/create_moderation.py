from openai import AsyncOpenAI
from utils.load_yaml import load_yaml

import json
import openai
import traceback

async def create_moderation(input_text):
    try:
        config = load_yaml(helpers.path_config_yaml)
        api_key = config['api_keys']['api_key_1']
        ai_client = AsyncOpenAI(api_key=api_key)
        response = await ai_client.moderations.create(
            model='omni-moderation-latest',
            input=input_text,
        )
        moderation = response.to_dict() if hasattr(response, 'to_dict') else response
        moderation_dict = json.loads(moderation)
        flagged = moderation_dict.get('results', [{}])[0].get('flagged', False)
        yield flagged
    except Exception as e:
        yield {'error': traceback.format_exc()}
