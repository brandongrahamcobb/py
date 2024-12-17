from utils.config import Config
from utils.increment_version import increment_version
from utils.setup_logging import setup_logging

import asyncio
import discord
import utils.helpers as helpers

async def main():
    config = Config().get_config()

    async with Discord(
        initial_extensions=config['discord_cogs'],
        intents=config['discord_intents'],
        testing_guild_id=config['discord_testing_guild_id'],
    ) as bot:
        bot.config = config
        increment_version(config, helpers.PATH_CONFIG_YAML)
        setup_logging(config, helpers.PATH_LOG)
        await bot.start(config['token'])

if __name__ == '__main__':
    asyncio.run(main())
