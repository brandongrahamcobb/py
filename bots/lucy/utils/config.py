''' config.py  The purpose of this program is to provide my primary configuaration script from cd ../.
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
from os import makedirs
from os.path import isfile, dirname
from typing import Dict, Any
from utils.load_yaml import load_yaml
from utils.prompt_for_values import prompt_for_values
from utils.setup_logging import logger

import utils.helpers as helpers
import yaml

class Config:
    _config = None  # Class variable to hold the config config

    @classmethod
    def get_config(cls) -> Dict[str, Any]:
        logger.info('Attempting to load configuration.')
        
        if cls._config is None:
            if isfile(helpers.PATH_CONFIG_YAML):
                logger.info(f'Config file found at {helpers.PATH_CONFIG_YAML}. Loading configuration.')
                config = load_yaml(helpers.PATH_CONFIG_YAML)
                
                # Ensure 'api_keys' exists and set defaults
                config['api_keys'] = config.get('api_keys', {})
                logger.info('API keys initialized.')

                # Loop through and prompt for each API key and associated values
                for i in range(1, 21):
                    key = f'api_key_{i}'
                    logger.info(f'Prompting for values for {key}.')
                    config['api_keys'][key] = config['api_keys'].get(key, {})
                    config['api_keys'][key]['api_key'] = prompt_for_values(f'Enter API key {i}', config['api_keys'][key].get('api_key', ''))
                    config['api_keys'][key]['client_id'] = prompt_for_values(f'Enter client ID for API key {i}', config['api_keys'][key].get('client_id', ''))
                    config['api_keys'][key]['client_secret'] = prompt_for_values(f'Enter client secret for API key {i}', config['api_keys'][key].get('client_secret', ''))
                    config['api_keys'][key]['redirect_uri'] = prompt_for_values(f'Enter redirect URI for API key {i}', config['api_keys'][key].get('redirect_uri', ''))

                # Prompt for additional configuration values
                logger.info('Prompting for other configuration values.')
                config['database_url'] = prompt_for_values('Enter the database URL', config.get('database_url', helpers.DATABASE_URL))
                config['discord_character_limit'] = prompt_for_values('Enter the character limit', config.get('discord_character_limit', helpers.DISCORD_CHARACTER_LIMIT))
                config['discord_cogs'] = prompt_for_values('Enter the cogs.', config.get('discord_cogs', helpers.DISCORD_COGS))
                config['discord_command_prefix'] = prompt_for_values('Enter the command prefix', config.get('discord_command_prefix', helpers.DISCORD_COMMAND_PREFIX))
                config['discord_intents'] = prompt_for_values('Enter the intents', config.get('discord_intents', helpers.DISCORD_INTENTS))
                config['discord_owner_id'] = prompt_for_values('Enter the owner ID', config.get('discord_owner_id', helpers.DISCORD_OWNER_ID))
                config['discord_role_pass'] = prompt_for_values('Vegan pass?', config.get('discord_role_pass', helpers.DISCORD_ROLE_PASS))
                config['discord_testing_guild_id'] = prompt_for_values('Enter the testing guild ID', config.get('discord_testing_guild_id', helpers.DISCORD_TESTING_GUILD_ID))
                config['discord_token'] = prompt_for_values('Enter the bot token', config.get('discord_token', ''))
                config['logging_level'] = prompt_for_values('Enter the logging level', config.get('logging_level', helpers.LOGGING_LEVEL))
                config['openai_chat_add_completion_to_history'] = prompt_for_values(f'Add completions to conversations?', config.get('openai_chat_add_completion_to_history', helpers.OPENAI_CHAT_ADD_COMPLETION_TO_HISTORY))
                config['openai_chat_max_tokens'] = prompt_for_values(f'Enter max tokens for completion', config.get('openai_chat_max_tokens', helpers.OPENAI_CHAT_MAX_TOKENS))
                config['openai_chat_moderation_model'] = prompt_for_values(f'Enter the chat model for moderation: {helpers.OPENAI_CHAT_MODELS}', config.get('openai_chat_moderation_model', helpers.OPENAI_CHAT_MODERATION_MODEL))
                config['openai_chat_model'] = prompt_for_values(f'Enter the chat model for ChatGPT: {helpers.OPENAI_CHAT_MODELS}', config.get('openai_chat_model', helpers.OPENAI_CHAT_MODEL))
                config['openai_chat_moderation'] = prompt_for_values(f'Enable or disable text moderation (True/False)', config.get('openai_moderation_text', helpers.OPENAI_CHAT_MODERATION))
                config['openai_chat_moderation_sys_input'] = prompt_for_values(f'Enter the text moderation system input', config.get('openai_chat_moderation_sys_input', helpers.OPENAI_CHAT_MODERATION_SYS_INPUT))
                config['openai_chat_n'] = prompt_for_values(f'Enter the chat response count', config.get('openai_n', helpers.OPENAI_CHAT_N))
                config['openai_chat_response_format'] = prompt_for_values(f'Enter the chat response format', config.get('openai_chat_response_format', helpers.OPENAI_CHAT_RESPONSE_FORMAT))
                config['openai_chat_store'] = prompt_for_values(f'Enable or disable chat storage (True/False)', config.get('openai_chat_store', helpers.OPENAI_CHAT_STORE))
                config['openai_chat_stream'] = prompt_for_values(f'Enable or disable chat streaming (True/False)', config.get('openai_chat_stream', helpers.OPENAI_CHAT_STREAM))
                config['openai_chat_stop'] = prompt_for_values(f'Enter the stop chat criteria', config.get('openai_chat_stop', helpers.OPENAI_CHAT_STOP))
                default_chat_sys_input = helpers.OPENAI_CHAT_SYS_INPUT if helpers.OPENAI_CHAT_MODEL in helpers.OPENAI_CHAT_MODELS['deprecated'] else ''
                config['openai_chat_sys_input'] = prompt_for_values('Enter the chat system input', config.get('openai_chat_sys_input', default_chat_sys_input))
                config['openai_chat_temperature'] = prompt_for_values(f'Enter the chat temperature from 0.0 to 2.0', config.get('openai_chat_temperature', helpers.OPENAI_CHAT_TEMPERATURE))
                config['openai_chat_top_p'] = prompt_for_values(f'Enter the chat top p', config.get('openai_chat_top_p', helpers.OPENAI_CHAT_TOP_P))
                config['openai_chat_use_history'] = prompt_for_values(f'Use chat history?', config.get('openai_chat_use_history', helpers.OPENAI_CHAT_USE_HISTORY))
                config['openai_chat_user'] = prompt_for_values(f'Enter your chat name', config.get('openai_chat_user', helpers.OPENAI_CHAT_USER))
                config['openai_moderation_image'] = prompt_for_values(f'Enable or disable image moderation (True/False)', config.get('openai_moderation_image', helpers.OPENAI_MODERATION_IMAGE))
                config['openai_moderation_model'] = prompt_for_values(f'Enter the image moderation model: {helpers.OPENAI_MODERATION_MODEL}', config.get('openai_moderation_model', helpers.OPENAI_MODERATION_MODEL))
                config['openai_moderation_warning'] = prompt_for_values(f'Enter the moderation warning', config.get('openai_moderation_warning', helpers.OPENAI_MODERATION_WARNING))
                config['openai_organization'] = prompt_for_values(f'Enter the organization ID', config.get('openai_organization', helpers.OPENAI_CHAT_HEADERS['OpenAI-Organization']))
                config['openai_project'] = prompt_for_values(f'Enter the project ID', config.get('openai_project', helpers.OPENAI_CHAT_HEADERS['OpenAI-Project']))
                config['twitch_token'] = prompt_for_values(f'Enter the Twitch token', config.get('twitch_token', ''))
                config['user_agent'] = prompt_for_values('Enter the User-Agent header', config.get('user_agent', helpers.USER_AGENT))
                config['version'] = prompt_for_values('Enter the bot version', config.get('version', helpers.VERSION))

                logger.info('Configuration loaded and all values prompted.')
            else:
                logger.info(f'Config file not found. Creating new config file at {helpers.PATH_CONFIG_YAML}.')
                makedirs(dirname(helpers.PATH_CONFIG_YAML), exist_ok=True)
                config = {
                    'api_keys': {
                        f'api_key_{i}': {
                            'api_key': prompt_for_values(f'Enter API key {i}', ''),
                            'client_id': prompt_for_values(f'Enter client ID for API key {i}', ''),
                            'client_secret': prompt_for_values(f'Enter client secret for API key {i}', ''),
                            'redirect_uri': prompt_for_values(f'Enter redirect URI for API key {i}', '')
                        } for i in range(1, 21)
                    },
                    'database_url': prompt_for_values('Enter the database URL', helpers.DATABASE_URL),
                    'discord_character_limit': prompt_for_values('Enter the character limit', helpers.DISCORD_CHARACTER_LIMIT),
                    'discord_cogs': prompt_for_values('Enter the cogs.', helpers.DISCORD_COGS),
                    'discord_command_prefix': prompt_for_values('Enter the command prefix', helpers.DISCORD_COMMAND_PREFIX),
                    'discord_intents': prompt_for_values('Enter the intents', helpers.DISCORD_INTENTS),
                    'discord_owner_id': prompt_for_values('Enter the owner ID', helpers.DISCORD_OWNER_ID),
                    'discord_role_pass': prompt_for_values('Vegan pass?', helpers.DISCORD_ROLE_PASS),
                    'discord_testing_guild_id': prompt_for_values('Enter the testing guild ID', helpers.DISCORD_TESTING_GUILD_ID),
                    'discord_token': prompt_for_values('Enter the bot token', ''),
                    'logging_level': prompt_for_values('Enter the logging level', helpers.LOGGING_LEVEL),
                    'openai_chat_add_completion_to_history': prompt_for_values(f'Add completions to conversations?', helpers.OPENAI_CHAT_ADD_COMPLETION_TO_HISTORY),
                    'openai_chat_max_tokens': prompt_for_values(f'Enter max tokens for completion', helpers.OPENAI_CHAT_MAX_TOKENS),
                    'openai_chat_model': prompt_for_values(f'Enter the chat model for ChatGPT: {helpers.OPENAI_CHAT_MODELS}', helpers.OPENAI_CHAT_MODEL),
                    'openai_chat_moderation_model': prompt_for_values(f'Enter the chat model for moderation: {helpers.OPENAI_CHAT_MODELS}', helpers.OPENAI_CHAT_MODERATION_MODEL),
                    'openai_chat_moderation': prompt_for_values('Enable or disable text moderation (True/False)', helpers.OPENAI_CHAT_MODERATION),
                    'openai_chat_n': prompt_for_values('Enter the response count', helpers.OPENAI_CHAT_N),
                    'openai_chat_response_format': prompt_for_values(f'Enter the chat response format', helpers.OPENAI_CHAT_RESPONSE_FORMAT),
                    'openai_chat_store': prompt_for_values('Enable or disable storage (True/False)', helpers.OPENAI_CHAT_STORE),
                    'openai_chat_stream': prompt_for_values('Enable or disable streaming (True/False)', helpers.OPENAI_CHAT_STREAM),
                    'openai_chat_stop': prompt_for_values('Enter the stop criteria', helpers.OPENAI_CHAT_STOP),
                    'openai_chat_sys_input': prompt_for_values('What is the system_input?', helpers.OPENAI_CHAT_SYS_INPUT) if helpers.OPENAI_CHAT_MODEL in helpers.OPENAI_CHAT_MODELS['deprecated'] else '',
                    'openai_chat_temperature': prompt_for_values('Enter the temperature (0.0 to 2.0)', helpers.OPENAI_CHAT_TEMPERATURE),
                    'openai_chat_top_p': prompt_for_values(f'Enter the chat top p', helpers.OPENAI_CHAT_TOP_P),
                    'openai_chat_add_completion_to_history': prompt_for_values(f'Add completions to conversations?', helpers.OPENAI_CHAT_ADD_COMPLETION_TO_HISTORY),
                    'openai_chat_use_history': prompt_for_values(f'Use chat history?', helpers.OPENAI_CHAT_USE_HISTORY),
                    'openai_chat_user': prompt_for_values('Enter your name', helpers.OPENAI_CHAT_USER),
                    'openai_moderation_image': prompt_for_values('Enable or disable image moderation (True/False)', helpers.OPENAI_MODERATION_IMAGE),
                    'openai_moderation_model': prompt_for_values(f'Enter the image moderation model: {helpers.OPENAI_MODERATION_MODEL}', helpers.OPENAI_MODERATION_MODEL),
                    'openai_moderation_warning': prompt_for_values('Enter the moderation warning', helpers.OPENAI_MODERATION_WARNING),
                    'openai_organization': prompt_for_values('Enter the organization ID', helpers.OPENAI_CHAT_HEADERS['OpenAI-Organization']),
                    'openai_project': prompt_for_values('Enter the project ID', helpers.OPENAI_CHAT_HEADERS['OpenAI-Project']),
                    'twitch_token': prompt_for_values(f'Enter the Twitch token', ''),
                    'user_agent': prompt_for_values('Enter the User-Agent header', helpers.USER_AGENT),
                    'version': prompt_for_values('Enter the bot version', helpers.VERSION)
                }

                with open(helpers.PATH_CONFIG_YAML, 'w') as file:
                    yaml.dump(config, file)
                    logger.info(f'Configuration saved to {helpers.PATH_CONFIG_YAML}.')

            cls._config = config
            logger.info('Configuration successfully loaded.')
        else:
            logger.info('Configuration already loaded, returning cached version.')
        
        return cls._config
