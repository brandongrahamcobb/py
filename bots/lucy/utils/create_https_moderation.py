''' create_https_moderation.py Still, OpenAI's python SDK is way better than this. However,
                               I'm proud of this program. It's my current work of art.
                               Functioning from cd ../
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
from utils.nlp_utils import NLPUtils

import aiohttp
import datetime
import json
import openai
import traceback
import utils.helpers as helpers

async def create_https_moderation(custom_id, input_array, model):
    try:
        config = load_yaml(helpers.PATH_CONFIG_YAML)
        api_key = config['api_keys']['api_key_1']
        ai_client = AsyncOpenAI(api_key=api_key)
        headers = {}
        headers.update({'Authorization': f'Bearer {api_key}'})
        request_data = {
            'input': input_array,
            'model': model,
        }
        async with aiohttp.ClientSession() as session:
            try:
                async with session.post(url=helpers.OPENAI_ENDPOINT_URLS['moderations'], headers=headers, json=request_data) as moderation_object:
                    if moderation_object.status == 200:
                        response_data = await moderation_object.json()
                        yield response_data
                    else:
                        yield {'error': await moderation_object.text()}
                    yield await moderation_object.json()
            except Exception as e:
                yield traceback.format_exc()
    except Exception as e:
        yield traceback.format_exc()

