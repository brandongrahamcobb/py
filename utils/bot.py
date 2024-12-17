from vyrtuous import Vyrtuous

import discord
import increment_version
import setup_logging

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
        bot.helpers = helpers
        increment_version(config)
        setup_logging(config)
        await bot.start(config['token'])
