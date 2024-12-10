from quart import Quart, request, Response, render_template
import asyncio
from openai import AsyncOpenAI
import yaml
import uuid
import logging
from collections import defaultdict
import os
import traceback

# Instantiate the Quart app
app = Quart(__name__)

sys_models = ['gpt-3.5-turbo', 'gpt-4', 'gpt-4-turbo', 'gpt-4o', 'gpt-4o-mini']
other_models = ['o1', 'o1-mini', 'o1-preview']

def prompt_for_values(prompt: str, default_value: str) -> str:
    """Prompt user for input, with a default value if not provided."""
    value = input(f'{prompt} [{default_value}]: ')
    return value if value else default_value

class Configuration:
    _config = None

    def __init__(self, api_key: str, model: str, stream: bool, **kwargs):
        """Initialize the AI instance with provided parameters."""
        self.api_key = api_key
        self.model = model
        self.stream = stream
        super().__init__(**kwargs)

    @classmethod
    def _get_config(cls) -> dict:
        """Load configuration from 'openai.yaml' or prompt user for missing values."""
        if cls._config is None:
            config_path = 'openai.yaml'
            if os.path.isfile(config_path):
                with open(config_path, 'r') as file:
                    data = yaml.safe_load(file)
            else:
                data = {}

            # Fetch configuration values from user or file
            data['api_key'] = prompt_for_values('Enter your OpenAI API key', data.get('api_key', ''))
            data['model'] = prompt_for_values('Enter model name', data.get('model', 'o1-preview'))
            while (data['model'] not in sys_models) and (data['model'] not in other_models):
                data['model'] = prompt_for_values(
                    f'Enter model name: {sys_models} or {other_models}',
                    data.get('model', 'o1-preview')
                )

            if data['model'] in sys_models:
                data['sys_input'] = prompt_for_values('Enter sys_input, if any.', data.get('sys_input', ''))
            else:
                data['sys_input'] = None
            data['stream'] = prompt_for_values(
                'Enable streaming? (True/False)',
                str(data.get('stream', True))
            ).lower() == 'true'

            data['logging_level'] = prompt_for_values('Enter the logging level', data.get('logging_level', 'INFO'))
            data['version'] = prompt_for_values('Enter the AI version', data.get('version', '1.0.0'))

            # Save to file
            with open(config_path, 'w') as file:
                yaml.dump(data, file)

            cls._config = data
        return cls._config

    @classmethod
    async def load_config(cls) -> dict:
        """Return the configuration."""
        return cls._get_config()

class GPT:
    def __init__(self, api_key: str, model: str, stream: bool):
        self.client = AsyncOpenAI(api_key=api_key)
        self.conversations = defaultdict(list)
        self.model = model
        self.stream = stream

    async def completion(self, input_text: str, conversation_id: str, sys_input: str = None):
        """Generate completion using OpenAI's API and stream the result."""
        try:
            messages = self.conversations[conversation_id]

            if sys_input and not any(msg['role'] == 'system' for msg in messages):
                messages.insert(0, {'role': 'system', 'content': sys_input})

            messages.append({'role': 'user', 'content': input_text})

            response_stream = await self.client.chat.completions.create(
                model=self.model,
                messages=messages,
                stream=self.stream
            )
            async def stream_response():
                full_response = ''
                async for chunk in response_stream:
                    delta = chunk.choices[0].delta
                    content = delta.get('content', '')
                    if content:
                        full_response += content
                        yield content
                # Append assistant's full response to the conversation history
                self.conversations[conversation_id].append({'role': 'assistant', 'content': full_response})

            headers = {
                'Content-Type': 'text/plain',
                'Cache-Control': 'no-cache',
                'X-Accel-Buffering': 'no'  # Disable buffering for Nginx
            }
            return Response(stream_response(), headers=headers)
        except Exception as e:
            logging.error(f"Error occurred: {traceback.format_exc()}")
            return Response(f"Error: {e}", status=500)

# Load configuration and initialize GPT instance
config = Configuration._get_config()
logging.basicConfig(level=config.get('logging_level', 'INFO'))
gpt = GPT(api_key=config['api_key'], model=config['model'], stream=config['stream'])

@app.route('/api/main', methods=['POST'])
async def api_main():
    data = await request.get_json()
    user_input = data.get('input', '')
    if not user_input:
        return {'error': 'No input provided'}, 400

    conversation_id = data.get('conversation_id') or str(uuid.uuid4())
    return await gpt.completion(user_input, conversation_id, sys_input=config['sys_input'])

@app.route('/')
async def index():
    return await render_template('VeganAI.html')

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
