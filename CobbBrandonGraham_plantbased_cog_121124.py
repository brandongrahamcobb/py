''' my_cog.py
    Copyright (C) 2024 github.com/<your-github>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''
from discord.ext import commands
from recipe_scrapers import scrape_me

import discord

class PlantBasedCog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot

    async def fetch_recipe(ctx, url: str):
        try:
            scraper = scrape_me(url)
            title = scraper.title()
            total_time = scraper.total_time()
            ingredients = scraper.ingredients()
            instructions = scraper.instructions()
            embed = discord.Embed(title=title, color=discord.Color.green())
            if total_time:
                embed.add_field(name="Total Time", value=f"{total_time} minutes", inline=False)
            embed.add_field(name="Ingredients", value='\n'.join(ingredients), inline=False)
            embed.add_field(name="Instructions", value=instructions, inline=False)
            await ctx.send(embed=embed)
        except Exception as e:
            await ctx.send(f"An error occurred while fetching the recipe: {e}")


async def setup(bot: commands.Bot):
    await bot.add_cog(PlantBasedCog(bot))
