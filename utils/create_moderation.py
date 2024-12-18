from openai import AsyncOpenAI
from utils.load_yaml import load_yaml

import openai
import traceback
import utils.helpers as helpers

async def create_moderation(input_text):
    try:
        config = load_yaml(helpers.PATH_CONFIG_YAML)
        api_key = config['api_keys']['api_key_1']
        ai_client = AsyncOpenAI(api_key=api_key)
        response = await ai_client.moderations.create(
            model='omni-moderation-latest',
            input=input_text,
        )
        yield await response.json()
    except Exception as e:
        yield {'error': traceback.format_exc()}
