from discord.ext import commands
from typing import List, Optional

import asyncio
import discord
import json
import logging
import logging.handlers
import os
import requests

class Lucy(commands.Bot):

    def __init__(
        self,
        *args,
        initial_extensions: List[str],
        testing_guild_id: Optional[int] = None,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.initial_extensions = initial_extensions
        self.testing_guild_id = testing_guild_id

    async def setup_hook(self) -> None:
        for extension in self.initial_extensions:
            await self.load_extension(extension)
        if self.testing_guild_id:
            guild = discord.Object(self.testing_guild_id)
            self.tree.copy_global_to(guild = guild)
            await self.tree.sync(guild = guild)
