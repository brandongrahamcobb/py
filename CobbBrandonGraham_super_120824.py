
#    Lucy, a discord.py bot, is an open-source package containing four cogs and a helper file.
#    Copyright (C) 2024  Cobb, Brandon Graham

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

from os import getenv, makedirs
from os.path import abspath, dirname, exists, expanduser, join

import datetime as dt
import json
import logging
import logging.handlers
import os

global logger

class Lucy(commands.Bot):

    def __init__(
        self,
        *args,
        **kwargs
    ):
        self.CURRENT_DATE = dt.datetime.now().strftime('%d%m%y')
        self.COBBBRANDONGRAHAM = expanduser('~')
        if os.name == 'nt':
            self.LICENSE = join(self.COBBBRANDONGRAHAM, 'Program Files', 'lucy-package', 'src', f'lucy_{self.CURRENT_DATE}', 'LICENSE') # ON PURPOSE 'LICENSE' NOT 'self.PATH_DEST_LICENSE'
            self.PATH_DEST_JSON= join(getenv('APPDATA'), '.config', 'CobbBrandonGraham', 'config.json')
            self.PATH_DEST_LOG = join(get(env('APPDATA'), '.log', 'discord.log')
            self.PATH_DEST_PY_ADMIN_COG = join(self.COBBBRANDONGRAHAM, 'Program Files', 'lucy-package', 'src', f'lucy_{self.CURRENT_DATE}', 'bot', 'cogs', 'admin_cog.py')
            self.PATH_DEST_PY_COGS_INIT = join(self.COBBBRANDONGRAHAM, 'Program Files', 'lucy-package', 'src', f'lucy_{self.CURRENT_DATE}', 'bot', 'cogs', '__init__.py')
            self.PATH_DEST_PY_GAME_COG = join(self.COBBBRANDONGRAHAM, 'Program Files', 'lucy-package', 'src', f'lucy_{self.CURRENT_DATE}', 'bot', 'cogs', 'game_cog.py')
            self.PATH_DEST_PY_HELPERS = join(self.COBBBRANDONGRAHAM, 'Program Files', 'lucy-package', 'src', f'lucy_{self.CURRENT_DATE}', 'bot', 'utils', 'helpers.py')
            self.PATH_DEST_PY_INIT = join(self.COBBBRANDONGRAHAM, 'Program Files', 'lucy-package', 'src', f'lucy_{self.CURRENT_DATE}', 'bot', '__init__.py')
            self.PATH_DEST_PY_MAIN = join(self.COBBBRANDONGRAHAM, 'Program Files', 'lucy-package', 'src', f'lucy_{self.CURRENT_DATE}', 'bot', 'main.py')
            self.PATH_DEST_PY_MY_COG = join(self.COBBBRANDONGRAHAM, 'Program Files', 'lucy-package', 'src', f'lucy_{self.CURRENT_DATE}', 'bot', 'cogs', 'my_cog.py')
            self.PATH_DEST_PY_USER_COG = join(self.COBBBRANDONGRAHAM, 'Program Files', 'lucy-package', 'src', f'lucy_{self.CURRENT_DATE}', 'bot', 'cogs', 'user_cog.py')
            self.PATH_DEST_PY_UTILS_INIT = join(self.COBBBRANDONGRAHAM, 'Program Files', 'lucy-package', 'src', f'lucy_{self.CURRENT_DATE}', 'bot', 'utils', '__init__.py')
            self.PATH_DEST_TXT_REQUIREMENTS = join(self.COBBBRANDONGRAHAM, 'Program Files', 'lucy-package', 'src', f'lucy_{self.CURRENT_DATE}', 'requirements.txt')
            self.README = join(self.COBBBRANDONGRAHAM, 'Program Files', 'lucy-package', 'src', f'lucy_{self.CURRENT_DATE}', 'README') # ON PURPOSE 'README' NOT 'self.PATH_DEST_README'
        else:
            self.LICENSE = join('usr', 'share', 'lucy',  'LICENSE') # ON PURPOSE 'LICENSE' NOT 'self.PATH_DEST_LICENSE'
            self.PATH_DEST_JSON = join(self.COBBBRANDONGRAHAM, '.config', 'CobbBrandonGraham', 'config.json')
            self.PATH_DEST_LOG = join(self.COBBBRANDONGRAHAM, '.log', 'discord.log')
            self.PATH_DEST_PY_ADMIN_COG = join('usr', 'share', 'lucy', 'bot', 'cogs', 'admin_cog.py')
            self.PATH_DEST_PY_COGS_INIT = join('usr', 'share', 'lucy', 'bot', 'cogs', '__init__.py')
            self.PATH_DEST_PY_GAME_COG = join('usr', 'share', 'lucy', 'bot', 'cogs', 'game_cog.py')
            self.PATH_DEST_PY_HELPERS = join('usr', 'share', 'lucy', 'bot', 'utils', 'helpers.py')
            self.PATH_DEST_PY_INIT = join('usr', 'share', 'lucy', 'bot', '__init__.py')
            self.PATH_DEST_PY_MAIN = join('usr', 'share', 'lucy', 'bot', 'main.py')
            self.PATH_DEST_PY_MY_COG = join('usr', 'share', 'lucy', 'bot', 'cogs', 'my_cog.py')
            self.PATH_DEST_PY_USER_COG = join('usr', 'share', 'lucy', 'bot', 'cogs', 'user_cog.py')
            self.PATH_DEST_PY_UTILS_INIT = join('usr', 'share', 'lucy', 'bot', 'utils', '__init__.py')
            self.PATH_DEST_TXT_REQUIREMENTS = join('usr', 'share', 'lucy', 'requirements.txt')
            self.README = join('usr', 'share', 'lucy',  'README') # ON PURPOSE 'README' NOT 'self.PATH_DEST_README'
        config = Lucy.update()
        self.initial_extensions = config['initial_extensions']
        self.testing_guild_id = config['testing_guild_id']

    def detect_cogs(self):
        return [f for f in os.listdir(dirname(PATH_DEST_PY_COGS_INIT)) if f.endswith('.py') and not f.startswith('__')]

    async def setup_hook(self) -> None:
        for cog in self.initial_extensions:
            await self.load_extention(cog)
        if self.testing_guild_id:
            guild = discord.Object(self.testing_guild_id)
            self.tree.copy_global_to(guild=guild)
            await self.tree.sync(guild=guild)

    def setup_logging(self, logging_level):
        logging_level = logging_level.upper()
        logging.basicConfig(level=getattr(logging, logging_level))
        if not exists(self.PATH_DEST_LOG):
            os.makedirs(self.PATH_DEST_LOG)
        if not exists(self.PATH_DEST_LOG):
            open(self.PATH_DEST_LOG, 'a').close()
        file_handler = logging.FileHandler(self.PATH_DEST_LOG)
        file_handler.setLevel(logging_level)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        logger = logging.getLogger()
        logger.setLevel(logging_level)
        logger.addHandler(file_handler)

   def update(self):
        if os.path.isfile(self.PATH_DEST_JSON):
            with open(self.PATH_DEST_JSON, 'r') as file:
                data = json.load(file)
            data['cogs'] = self.prompt_for_values('Enter the cogs.', data.get('cogs', [
                'bot.cogs.admin_cog',
                'bot.cogs.game_cog',
                'bot.cogs.my_cog',
                'bot.cogs.user_cog'
            ]))
            data['command_prefix'] = self.prompt_for_values('Enter the command prefix', data.get('command_prefix', '!'))
            data['database_url'] = self.prompt_for_values('Enter the database URL', data.get('database_url', 'brandongcobb.com'))
            data['intents'] = self.prompt_for_values('Enter the intents', data.get('intents', 'discord.Intents.all()'))
            data['logging_level'] = self.prompt_for_values('Enter the logging level', data.get('logging_level', 'INFO'))
            data['owner_id'] = self.prompt_for_values('Enter the owner ID', data.get('owner_id', '154749533429956608'))
            data['testing_guild_id'] = self.prompt_for_values('Enter the testing guild ID', data.get('testing_guild_id', '1217326055111655507'))
            data['token'] = self.prompt_for_values('Enter the bot token', data.get('token', ''))
            data['user_agent'] = self.prompt_for_values('Enter the User-Agent header', data.get('user_agent', 'Lucy'))
            data['version'] = self.prompt_for_values('Enter the bot version', data.get('version', '1.0.0'))
            data['xp_rate'] = self.prompt_for_values('Enter the XP rate', data.get('xp_rate', '1'))
            for i in range(1, 21):
                key = f'api_key_{i}'
                current_key = data.get('api_keys', {}).get(key, '')
                data['api_keys'][key] = self.prompt_for_values(f'Enter API key {i}', current_key)
        else:
            data = {
                'api_keys': {},
                'cogs': self.prompt_for_values('Enter the cogs.', ['bot.cogs.admin_cog', 'bot.cogs.game_cog,', 'bot.cogs.my_cog', 'bot.cogs.user_cog'])
                'command_prefix': self.prompt_for_values('Enter the command prefix', '!'),
                'database_url': self.prompt_for_values('Enter the database URL', 'brandongcobb.com'),
                'intents': self.prompt_for_values('Enter the intents', 'discord.Intents.all()'),
                'logging_level': self.prompt_for_values('Enter the logging level', 'INFO'),
                'owner_id': self.prompt_for_values('Enter the your user ID', '154749533429956608'),
                'testing_guild_id': self.prompt_for_values('Enter the testing guild ID', '1217326055111655507'),
                'token': self.prompt_for_values('Enter the bot token', ''),
                'version': self.prompt_for_values('Enter the bot version', '1.0.0'),
                'user_agent': self.prompt_for_values('Enter the User-Agent header', 'Lucy'),
                'xp_rate': self.prompt_for_values('Enter the XP rate', '1'),
            ]'))
            }
            for i in range(1, 21):
                data['api_keys'][f'api_key_{i}'] = self.prompt_for_values(f'Enter API key {i}', '')
        current_version = data['version']
        major, minor, patch = map(int, current_version.split('.'))
        patch += 1
        if patch >= 10:
            patch = 0
            minor += 1
        if minor >= 10:
            minor = 0
            major += 1
        new_version = f'{major}.{minor}.{patch}'
        data['version'] = new_version
        with open(self.PATH_DEST_JSON, 'w') as file:
            json.dump(data, file, indent=4)
        return data

    def prompt_for_values(self, prompt, default_value):
        value = input(f'{prompt} [{default_value}]: ')
        return value if value else default_value

    def run_makepkg():
        try:
            result = subprocess.run(['makepkg', '-si'], check=True, text=True, capture_output=True)
            print('Command Output:', result.stdout)
            print('Command Error:', result.stderr)
        except subprocess.CalledProcessError as e:
            print('An error occurred while running makepkg:')
            print('Error Code:', e.returncode)
            print('Error Output:', e.stderr)

if __name__ is '__main__':

    lucy = Lucy()

    lucy.run_makepkg(self)

    config = lucy.update(self):

    lucy.setup_logging(self, config['logging_level'])
    
    async with Lucy() as bot:
        await bot.start(config['token'])
