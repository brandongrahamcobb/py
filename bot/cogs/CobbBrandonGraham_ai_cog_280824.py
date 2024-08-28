from discord.ext import commands
from bot.utils.helpers import check_for_updates

class AICog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
#        self.bot.loop.create_task(check_for_updates(bot))

async def setup(bot):
    await bot.add_cog(AICog(bot))
