from utils.config import Config
from utils.discord import Discord
from utils.increment_version import increment_version
from utils.setup_logging import setup_logging

import asyncio
import asyncpg
import discord
import utils.helpers as helpers

async def main():
    config = Config().get_config()

    async with asyncpg.create_pool(database='lucy', user='postgres', command_timeout=30) as pool:
        async with Discord(
            command_prefix=config['discord_command_prefix'],
            db_pool=pool,
            initial_extensions=config['discord_cogs'],
            intents=eval(config['discord_intents']),
            testing_guild_id=config['discord_testing_guild_id'],
        ) as bot:
            bot.config = config
            increment_version(config, helpers.PATH_CONFIG_YAML)
            setup_logging(config, helpers.PATH_LOG)
            await bot.start(config['discord_token'])

if __name__ == '__main__':
    asyncio.run(main())
