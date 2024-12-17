from utils.increment_version import increment_version
from utils.setup_logging import setup_logging
from utils.vyrtuous import Vyrtuous

import asyncio
import discord
import utils.helpers as helpers

async def main():
    intents = discord.Intents.all()

    async with Vyrtuous(
        command_prefix='!',
        intents=intents,
        initial_extensions=[
            'bot.cogs.hybrid',
            'bot.cogs.indica',
            'bot.cogs.sativa'
        ],
    ) as bot:

        config = await bot.load_config()
        bot.config = config
        increment_version(config, helpers.PATH_CONFIG_YAML)
        setup_logging(config, helpers.PATH_LOG)
        await bot.start(config['token'])

if __name__ == '__main__':
    asyncio.run(main())
