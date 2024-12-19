from discord.ext import commands
from typing import List, Optional

import asyncpg
import discord

class Discord(commands.Bot):

    def __init__(
        self,
        *args,
        initial_extensions: List[str],
        db_pool: asyncpg.Pool,
        testing_guild_id: Optional[int] = None,
        **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.db_pool = db_pool
        self.testing_guild_id = testing_guild_id
        self.initial_extensions = initial_extensions

    async def setup_hook(self) -> None:
        for cog in self.initial_extensions:
            await self.load_extension(cog)
        if self.testing_guild_id:
            guild = discord.Object(self.testing_guild_id)
            self.tree.copy_global_to(guild=guild)
            await self.tree.sync(guild=guild)
