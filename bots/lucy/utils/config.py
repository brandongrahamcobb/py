from os import makedirs
from os.path import isfile, dirname
from typing import Dict, Any
import utils.helpers as helpers
from utils.load_yaml import load_yaml
from utils.prompt_for_values import prompt_for_values

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
                config['database_url'] = prompt_for_values('Enter the database URL', config.get('database_url', helpers.DATABASE_URL))
                config['discord_character_limit'] = prompt_for_values('Enter the character limit', config.get('discord_character_limit', helpers.DISCORD_CHARACTER_LIMIT))
                config['discord_cogs'] = prompt_for_values('Enter the cogs.', config.get('discord_cogs', helpers.DISCORD_COGS))
                config['discord_command_prefix'] = prompt_for_values('Enter the command prefix', config.get('discord_command_prefix', helpers.DISCORD_COMMAND_PREFIX))
                config['discord_intents'] = prompt_for_values('Enter the intents', config.get('discord_intents', helpers.DISCORD_INTENTS))
                config['discord_owner_id'] = prompt_for_values('Enter the owner ID', config.get('discord_owner_id', helpers.DISCORD_OWNER_ID))
                config['discord_testing_guild_id'] = prompt_for_values('Enter the testing guild ID', config.get('discord_testing_guild_id', helpers.DISCORD_TESTING_GUILD_ID))
                config['discord_token'] = prompt_for_values('Enter the bot token', config.get('discord_token', ''))
                config['logging_level'] = prompt_for_values('Enter the logging level', config.get('logging_level', helpers.LOGGING_LEVEL))
                config['openai_chat_max_tokens'] = prompt_for_values(f'Enter max tokens for completion', config.get('openai_chat_max_tokens', helpers.OPENAI_CHAT_MAX_TOKENS))
                config['openai_chat_moderation_model'] = prompt_for_values(f'Enter the chat model for moderation: {helpers.OPENAI_CHAT_MODELS}', config.get('openai_chat_moderation_model', helpers.OPENAI_CHAT_MODERATION_MODEL))
                config['openai_chat_model'] = prompt_for_values(f'Enter the chat model for ChatGPT: {helpers.OPENAI_CHAT_MODELS}', config.get('openai_chat_model', helpers.OPENAI_CHAT_MODEL))
                config['openai_chat_moderation'] = prompt_for_values(f'Enable or disable text moderation (True/False)', config.get('openai_moderation_text', helpers.OPENAI_CHAT_MODERATION))
                config['openai_chat_moderation_sys_input'] = prompt_for_values(f'Enter the text moderation system input', config.get('openai_chat_moderation_sys_input', helpers.OPENAI_CHAT_MODERATION_SYS_INPUT))
                config['openai_chat_n'] = prompt_for_values(f'Enter the chat response count', config.get('openai_n', helpers.OPENAI_CHAT_N))
                config['openai_fine_tuning_response_format'] = prompt_for_values(f'Enter the fine_tuning response format', config.get('openai_fine_tuning_response_format', helpers.OPENAI_FINE_TUNING_RESPONSE_FORMAT))
                config['openai_chat_response_format'] = prompt_for_values(f'Enter the chat response format', config.get('openai_chat_response_format', helpers.OPENAI_CHAT_RESPONSE_FORMAT))
                config['openai_chat_store'] = prompt_for_values(f'Enable or disable chat storage (True/False)', config.get('openai_chat_store', helpers.OPENAI_CHAT_STORE))
                config['openai_chat_stream'] = prompt_for_values(f'Enable or disable chat streaming (True/False)', config.get('openai_chat_stream', helpers.OPENAI_CHAT_STREAM))
                config['openai_chat_stop'] = prompt_for_values(f'Enter the stop chat criteria', config.get('openai_chat_stop', helpers.OPENAI_CHAT_STOP))
                default_chat_sys_input = helpers.OPENAI_CHAT_SYS_INPUT if helpers.OPENAI_CHAT_MODEL in helpers.OPENAI_CHAT_MODELS['deprecated'] else ''
                config['openai_chat_sys_input'] = prompt_for_values('Enter the chat system input', config.get('openai_chat_sys_input', default_chat_sys_input))
                config['openai_chat_temperature'] = prompt_for_values(f'Enter the chat temperature from 0.0 to 2.0', config.get('openai_chat_temperature', helpers.OPENAI_CHAT_TEMPERATURE))
                config['openai_chat_top_p'] = prompt_for_values(f'Enter the chat top p', config.get('openai_chat_top_p', helpers.OPENAI_CHAT_TOP_P))
                config['openai_chat_user'] = prompt_for_values(f'Enter your chat name', config.get('openai_chat_user', helpers.OPENAI_CHAT_USER))
                config['openai_moderation_image'] = prompt_for_values(f'Enable or disable image moderation (True/False)', config.get('openai_moderation_image', helpers.OPENAI_MODERATION_IMAGE))
                config['openai_moderation_model'] = prompt_for_values(f'Enter the image moderation model: {helpers.OPENAI_MODERATION_MODEL}', config.get('openai_moderation_image', helpers.OPENAI_MODERATION_MODEL))
                config['openai_moderation_warning'] = prompt_for_values(f'Enter the moderation warning', config.get('openai_moderation_warning', helpers.OPENAI_MODERATION_WARNING))
                config['openai_organization'] = prompt_for_values(f'Enter the organization ID', config.get('openai_organization', helpers.OPENAI_CHAT_HEADERS['OpenAI-Organization']))
                config['openai_project'] = prompt_for_values(f'Enter the project ID', config.get('openai_project', helpers.OPENAI_CHAT_HEADERS['OpenAI-Project']))
                config['user_agent'] = prompt_for_values('Enter the User-Agent header', config.get('user_agent', helpers.USER_AGENT))
                config['version'] = prompt_for_values('Enter the bot version', config.get('version', helpers.VERSION))
            else:
                makedirs(dirname(helpers.PATH_CONFIG_YAML), exist_ok=True)
                config = {
                    'api_keys': {f'api_key_{i}': prompt_for_values(f'Enter API key {i}', '') for i in range(1, 21)},
                    'database_url': prompt_for_values('Enter the database URL', helpers.DATABASE_URL),
                    'discord_character_limit': prompt_for_values('Enter the character limit', helpers.DISCORD_CHARACTER_LIMIT),
                    'discord_cogs': prompt_for_values('Enter the cogs.', helpers.DISCORD_COGS),
                    'discord_command_prefix': prompt_for_values('Enter the command prefix', helpers.DISCORD_COMMAND_PREFIX),
                    'discord_intents': prompt_for_values('Enter the intents', helpers.DISCORD_INTENTS),
                    'discord_owner_id': prompt_for_values('Enter the owner ID', helpers.DISCORD_OWNER_ID),
                    'discord_testing_guild_id': prompt_for_values('Enter the testing guild ID', helpers.DISCORD_TESTING_GUILD_ID),
                    'discord_token': prompt_for_values('Enter the bot token', ''),
                    'logging_level': prompt_for_values('Enter the logging level', helpers.LOGGING_LEVEL),
                    'openai_chat_max_tokens': prompt_for_values(f'Enter max tokens for completion', helpers.OPENAI_CHAT_MAX_TOKENS),
                    'openai_chat_model': prompt_for_values(f'Enter the chat model for ChatGPT: {helpers.OPENAI_CHAT_MODELS}', helpers.OPENAI_CHAT_MODEL),
                    'openai_chat_moderation_model': prompt_for_values(f'Enter the chat model for moderation: {helpers.OPENAI_CHAT_MODELS}', helpers.OPENAI_CHAT_MODERATION_MODEL),
                    'openai_chat_moderation': prompt_for_values('Enable or disable text moderation (True/False)', helpers.OPENAI_CHAT_MODERATION),
                    'openai_chat_n': prompt_for_values('Enter the response count', helpers.OPENAI_CHAT_N),
                    'openai_fine_tuning_response_format': prompt_for_values(f'Enter the fine_tuning response format', helpers.OPENAI_FINE_TUNING_RESPONSE_FORMAT),
                    'openai_chat_response_format': prompt_for_values(f'Enter the chat response format', helpers.OPENAI_CHAT_RESPONSE_FORMAT),
                    'openai_chat_store': prompt_for_values('Enable or disable storage (True/False)', helpers.OPENAI_CHAT_STORE),
                    'openai_chat_stream': prompt_for_values('Enable or disable streaming (True/False)', helpers.OPENAI_CHAT_STREAM),
                    'openai_chat_stop': prompt_for_values('Enter the stop criteria', helpers.OPENAI_CHAT_STOP),
                    'openai_chat_sys_input': prompt_for_values('What is the system_input?', helpers.OPENAI_CHAT_SYS_INPUT) if helpers.OPENAI_CHAT_MODEL in helpers.OPENAI_CHAT_MODELS['deprecated'] else '',
                    'openai_chat_temperature': prompt_for_values('Enter the temperature (0.0 to 2.0)', helpers.OPENAI_CHAT_TEMPERATURE),
                    'openai_chat_top_p': prompt_for_values(f'Enter the chat top p', helpers.OPENAI_CHAT_TOP_P),
                    'openai_chat_user': prompt_for_values('Enter your name', helpers.OPENAI_CHAT_USER),
                    'openai_moderation_image': prompt_for_values('Enable or disable image moderation (True/False)', helpers.OPENAI_MODERATION_IMAGE),
                    'openai_moderation_model': prompt_for_values(f'Enter the image moderation model: {helpers.OPENAI_MODERATION_MODEL}', helpers.OPENAI_MODERATION_MODEL),
                    'openai_moderation_warning': prompt_for_values('Enter the moderation warning', helpers.OPENAI_MODERATION_WARNING),
                    'openai_organization': prompt_for_values('Enter the organization ID', helpers.OPENAI_CHAT_HEADERS['OpenAI-Organization']),
                    'openai_project': prompt_for_values('Enter the project ID', helpers.OPENAI_CHAT_HEADERS['OpenAI-Project']),
                    'user_agent': prompt_for_values('Enter the User-Agent header', helpers.USER_AGENT),
                    'version': prompt_for_values('Enter the bot version', helpers.VERSION)
                }
            with open(helpers.PATH_CONFIG_YAML, 'w') as file:
                yaml.dump(config, file)
            cls._config = config
        return cls._config

