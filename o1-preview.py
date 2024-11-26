""" openai.py The purpose of this program is to entend OpenAI's reach to other platforms
    Copyright (C) 2024  https://github.com/brandongrahamcobb

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
"""
from collections import defaultdict
from openai import AsyncOpenAI
import os
import threading
import sys
import yaml
import logging
from typing import List, Optional, Dict, Any
from os import getenv, makedirs
from os.path import abspath, dirname, expanduser, exists, isfile, join
import traceback

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def clear_screen():
    """Clear the terminal screen."""
    if sys.platform == "win32":
        os.system('cls')  # For Windows
    else:
        os.system('clear')  # For Linux/macOS

def prompt_for_values(prompt: str, default_value: str) -> str:
    """ Prompt user for input, with a default value if not provided. """
    value = input(f'{prompt} [{default_value}]: ')
    return value if value else default_value

class GPT:
    def __init__(self, api_key: str, model: str, stream: bool):
        self.client = AsyncOpenAI(api_key=api_key)
        self.conversations = defaultdict(list)
        self.model = model
        self.stream = stream

    async def completion(self, input_text: str, conversation_id: str):
        """ Generate completion using OpenAI's API and stream the result. """
        try:
            messages = self.conversations[conversation_id]
            messages.append({'role': 'user', 'content': input_text})

            stream = await self.client.chat.completions.create(
                model=self.model,
                messages=messages,
                stream=self.stream
            )
            full_response = ''
            async for chunk in stream:
                content = chunk.choices[0].delta.content
                if content is not None:
                    full_response += content

            self.conversations[conversation_id].append({'role': 'assistant', 'content': full_response})
            yield full_response
        except Exception as e:
            logger.error(f"Error occurred: {traceback.format_exc()}")
            yield traceback.format_exc()

class AI:
    _config = None

    def __init__(self, api_key: str, model: str, stream: bool, **kwargs):
        """ Initialize the AI instance with provided parameters. """
        self.api_key = api_key
        self.model = model
        self.stream = stream
        super().__init__(**kwargs)

    @classmethod
    def _get_config(cls) -> Dict[str, Any]:
        """ Load configuration from 'openai.yaml' or prompt user for missing values. """
        if cls._config is None:
            config_path = 'openai.yaml'
            if isfile(config_path):
                with open(config_path, 'r') as file:
                    data = yaml.safe_load(file)
            else:
                data = {}

            # Fetch configuration values from user or file
            data['api_key'] = prompt_for_values('Enter your OpenAI API key', data.get('api_key', ''))
            data['model'] = prompt_for_values('Enter model name', data.get('model', 'gpt-o1-preview'))
            data['stream'] = prompt_for_values('Enable streaming? (True/False)', str(data.get('stream', True))).lower() == 'true'
            data['logging_level'] = prompt_for_values('Enter the logging level', data.get('logging_level', 'INFO'))
            data['version'] = prompt_for_values('Enter the AI version', data.get('version', '1.0.0'))

            # Save to file
            with open(config_path, 'w') as file:
                yaml.dump(data, file)

            cls._config = data

        return cls._config

    async def load_config(self) -> Dict[str, Any]:
        """ Return the configuration. """
        return self._get_config()

async def main():
    """ Main function that loads config and runs the AI process. """
    config = AI._get_config()

    # Set logging level based on configuration
    logging_level = getattr(logging, config['logging_level'].upper(), logging.INFO)
    logger.setLevel(logging_level)

    # Initialize GPT instance with the loaded configuration
    gpt_instance = GPT(
        api_key=config['api_key'],
        model=config['model'],
        stream=config['stream']
    )

    # Example completion request
    conversation_id = '1234'  # Dummy conversation ID
    async for response in gpt_instance.completion('Hello, AI!', conversation_id):
        print(response)
    
    while True:
        input_text = input('You: ')
        clear_screen()
        async for response in gpt_instance.completion(input_text, conversation_id):
            print(f'AI: {response}')

if __name__ == '__main__':
    # Run the main function using asyncio
    import asyncio
    asyncio.run(main())
