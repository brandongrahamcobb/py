from os import makedirs
from os.path import isfile, dirname
from typing import Dict, Any

from prompt_for_values import prompt_for_values
import helpers as helpers
from load_yaml import load_yaml
import yaml

class Config:
    _config = None  # Class variable to hold the config config

    @classmethod
    def get_config(cls) -> Dict[str, Any]:
        if cls._config is None:
            if isfile(helpers.PATH_CONFIG_YAML):
                config = load_yaml(helpers.PATH_CONFIG_YAML)
                config['api_keys'] = config.get('api_keys', {})
                for i in range(1, 21):
                    key = f'api_key_{i}'
                    current_key = config['api_keys'].get(key, '')
                    config['api_keys'][key] = prompt_for_values(f'Enter API key {i}', current_key)
                config['database_url'] = prompt_for_values('Enter the database URL', config.get('database_url', ''))
                config['discord_cogs'] = prompt_for_values('Enter the cogs.', config.get('cogs', [
                    'bot.cogs.hybrid',
                    'bot.cogs.indica',
                    'bot.cogs.sativa'
                ]))
                config['discord_command_prefix'] = prompt_for_values('Enter the command prefix', config.get('discord_command_prefix', '!'))
                config['discord_intents'] = prompt_for_values('Enter the intents', config.get('discord_intents', 'discord.Intents.all()'))
                config['discord_owner_id'] = prompt_for_values('Enter the owner ID', config.get('discord_owner_id', '154749533429956608'))
                config['discord_testing_guild_id'] = prompt_for_values('Enter the testing guild ID', config.get('discord_testing_guild_id', '1286307083976835103'))
                config['discord_token'] = prompt_for_values('Enter the bot token', config.get('discord_token', ''))
                config['logging_level'] = prompt_for_values('Enter the logging level', config.get('logging_level', 'INFO'))
                config['openai_model'] = prompt_for_value(f'Enter the model: {helpers.OPENAI_MODELS}', config.get('openai_model', 'gpt-4o-mini')
                config['openai_n'] = prompt_for_value(f'Enter the response number', config.get('openai_n', '1')
                config['openai_organization'] = prompt_for_value(f'Enter the organization ID', config.get('openai_organization', OPENAI_HEADERS['OpenAI-Organization'])
                config['openai_project'] = prompt_for_value(f'Enter the project ID', config.get('openai_project', OPENAI_HEADERS['OpenAI-Project'])
                config['openai_store'] = prompt_for_value(f'Enter the storage preference (True/False)?', config.get('openai_store', 'True')
                config['openai_stream'] = prompt_for_value(f'Enter the stream preference (True/False)?', config.get('openai_stream', 'True')
                config['openai_stop'] = prompt_for_value(f'Enter the stop criteria?', config.get('openai_stop', '')
                config['openai_temperature'] = prompt_for_value(f'Enter from 0.0 to 2.0 tje temperature?', config.get('openai_temperature', '1.0')
                config['openai_user'] = prompt_for_value(f'Enter your name?', config.get('openai_user', 'Brandon Graham Cobb')
                config['user_agent'] = prompt_for_values('Enter the User-Agent header', config.get('user_agent', 'Vyrtuous'))
                config['version'] = prompt_for_values('Enter the bot version', config.get('version', '1.0.0'))
                if config['openai_model'] in helpers.OPENAI_MODELS['deprecated']:
                    config['openai_system_input'] = prompt_for_value(f'What is the system_input?', config.get('openai_system_input', '')
            else:
                makedirs(dirname(helpers.PATH_CONFIG_YAML), exist_ok=True)
                config = {
                    'api_keys': {},
                    'database_url': prompt_for_values('Enter the database URL', ''),
                    'discord_cogs': prompt_for_values('Enter the cogs.', [
                        'bot.cogs.hybrid',
                        'bot.cogs.indica',
                        'bot.cogs.sativa'
                    ]),
                    'discord_command_prefix': prompt_for_values('Enter the command prefix', '!'),
                    'discord_intents': prompt_for_values('Enter the intents', 'discord.Intents.all()'),
                    'discord_owner_id': prompt_for_values('Enter the owner ID', '154749533429956608'),
                    'discord_testing_guild_id': prompt_for_values('Enter the testing guild ID', '1286307083976835103'),
                    'discord_token': prompt_for_values('Enter the bot token', ''),
                    'logging_level': prompt_for_values('Enter the logging level', 'INFO'),
                    'openai_model': prompt_for_values(f'Enter the model', 'gpt-4o-mini'),
                    'openai_n': prompt_for_values(f'Enter the response number', '1'),
                    'openai_organization': prompt_for_values(f'Enter the organization ID', helpers.OPENAI_HEADERS['OpenAI-Organization']),
                    'openai_project': prompt_for_values('Enter the project ID', helpers.OPENAI_HEADERS['OpenAI-Project']),
                    'openai_store': prompt_for_values(f'Enter the storage preference (True/False)', 'True'),
                    'openai_stream': prompt_for_values(f'Enter the stream preference (True/False)', 'True'),
                    'openai_stop': prompt_for_values(f'Enter the stop criteria?', ''),
                    'openai_temperature': prompt_for_values(f'Enter the temperature (0.0 to 2.0)', '1.0'),
                    'openai_user': prompt_for_values(f'Enter your name', 'Brandon Graham Cobb'),
                    'user_agent': prompt_for_values('Enter the User-Agent header', 'Vyrtuous'),
                    'version': prompt_for_values('Enter the bot version', '1.0.0'),
                    'openai_system_input': prompt_for_values('What is the system_input?', ''),
                }
                for i in range(1, 21):
                    config['api_keys'][f'api_key_{i}'] = prompt_for_values(f'Enter API key {i}', '')
            with open(helpers.PATH_CONFIG_YAML, 'w') as file:
                yaml.dump(config, file)
            cls._config = config
        return cls._config
=
