from openai import AsyncOpenAI
from utils.load_yaml import load_yaml
from utils.setup_logging import logger

import openai
import traceback
import utils.helpers as helpers

async def create_completion(input_array):
    try:
        config = load_yaml(helpers.PATH_CONFIG_YAML)
        api_key = config['api_keys']['api_key_1']['api_key']
        ai_client = AsyncOpenAI(api_key=api_key)
        response = await ai_client.chat.completions.create(
            model='gpt-4o-mini',
            messages=input_array,
            response_format=helpers.OPENAI_CHAT_RESPONSE_FORMAT
        )
        yield response.choices[0].message.content
    except Exception as e:
        yield {'error': traceback.format_exc()}
