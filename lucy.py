from bot import Lucy
from discord.ext import commands
from typing import List, Optional

import asyncio
import discord
import json
import logging
import logging.handlers
import os

class Main:

    async def main():
    
        logger = logging.getLogger('discord')
        logger.setLevel(logging.INFO)
        handler = logging.handlers.RotatingFileHandler(
            filename=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'log', 'discord.log'),
            encoding='utf-8',
            maxBytes=32 * 1024 * 1024,
             backupCount=5,
        )
        dt_fmt = '%Y-%m-%d %H:%M:%S'
        formatter = logging.Formatter('[{asctime}] [{levelname:<8} {name}: {message}', dt_fmt, style='{')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        exts = ['vegan_university']
        intents = discord.Intents.all()
        intents.message_content = True
    
        async with Lucy(
            commands.when_mentioned_or('!'),
            initial_extensions = exts,
            intents = intents,
        ) as bot:
            dir = os.path.dirname(os.path.abspath(__file__))
            updir = os.path.join(dir, '..')
            dir_json = os.path.join(updir, 'json')
            file_json = os.path.join(dir_json, 'config.json')
            with open(file_json, 'r') as f:
                data = json.load(f)
                token = data['token']
                await bot.start(token)

asyncio.run(Main.main())
