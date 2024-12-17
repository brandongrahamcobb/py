from discord import commands
from typing import Any, Dict, List, Optional

import discord

class Vyrtuous(commands.Bot):

    _config = None                         # This is to establish a class variable.

    def __init__(
        self,
        *args,
        initial_extensions: List[str],
        testing_guild_id: Optional[int] = None,
        **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.testing_guild_id = testing_guild_id
        self.initial_extensions = initial_extensions

    @classmethod
    def _get_config(cls) -> Dict[str, Any]:
        if cls._config is None:
            if isfile(helpers.path_config_yaml):
                data = load_yaml(helpers.path_config_yaml)
                data['cogs'] = prompt_for_values('Enter the cogs.', data.get('cogs', [
                    'bot.cogs.hybrid',
                    'bot.cogs.indica',
                    'bot.cogs.sativa'
                ]))
                data['command_prefix'] = prompt_for_values('Enter the command prefix', data.get('command_prefix', '!'))
                data['database_url'] = prompt_for_values('Enter the database URL', data.get('database_url', 'brandongcobb.com'))
                data['intents'] = prompt_for_values('Enter the intents', data.get('intents', 'discord.Intents.all()'))
                data['logging_level'] = prompt_for_values('Enter the logging level', data.get('logging_level', 'INFO'))
                data['owner_id'] = prompt_for_values('Enter the owner ID', data.get('owner_id', '154749533429956608'))
                data['testing_guild_id'] = prompt_for_values('Enter the testing guild ID', data.get('testing_guild_id', '1286307083976835103'))
                data['token'] = prompt_for_values('Enter the bot token', data.get('token', ''))
                data['user_agent'] = prompt_for_values('Enter the User-Agent header', data.get('user_agent', 'Vyrtuous'))
                data['version'] = prompt_for_values('Enter the bot version', data.get('version', '1.0.0'))
                data['api_keys'] = data.get('api_keys', {})
                for i in range(1, 21):
                    key = f'api_key_{i}'
                    current_key = data['api_keys'].get(key, '')
                    data['api_keys'][key] = prompt_for_values(f'Enter API key {i}', current_key)
            else:
                makedirs(dirname(helpers.path_config_yaml), exist_ok=True)
                data = {
                    'api_keys': {},
                    'cogs': prompt_for_values('Enter the cogs.', [
                        'bot.cogs.hybrid',
                        'bot.cogs.indica',
                        'bot.cogs.sativa'
                    ]),
                    'command_prefix': prompt_for_values('Enter the command prefix', '!'),
                    'database_url': prompt_for_values('Enter the database URL', 'brandongcobb.com'),
                    'intents': prompt_for_values('Enter the intents', 'discord.Intents.all()'),
                    'logging_level': prompt_for_values('Enter the logging level', 'INFO'),
                    'owner_id': prompt_for_values('Enter the your user ID', '154749533429956608'),
                    'testing_guild_id': prompt_for_values('Enter the testing guild ID', '1286307083976835103'),
                    'token': prompt_for_values('Enter the bot token', ''),
                    'version': prompt_for_values('Enter the bot version', '1.0.0'),
                    'user_agent': prompt_for_values('Enter the User-Agent header', 'Vyrtuous'),
                }
                for i in range(1, 21):
                    data['api_keys'][f'api_key_{i}'] = prompt_for_values(f'Enter API key {i}', '')
            with open(helpers.path_config_yaml, 'w') as file:
                yaml.dump(data, file)
            cls._config = data
            return cls._config

    async def load_config(self) -> Dict[str, Any]:
        return self._get_config()

    async def setup_hook(self) -> None:
        for cog in self.initial_extensions:
            await self.load_extension(cog)
        if self.testing_guild_id:
            guild = discord.Object(self.testing_guild_id)
            self.tree.copy_global_to(guild=guild)
            await self.tree.sync(guild=guild)
