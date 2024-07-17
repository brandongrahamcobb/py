# cogs/reload.py

import discord
import os
import re

from discord.ext import commands

class Reload(commands.Cog):

    def __init__(self, bot):
        self.bot = bot
        self._last_member = None

    @commands.group(name='reload', hidden=True, invoke_without_command=True)
    async def _reload(self, ctx: commands.Context, *, module: str):
        """Reloads a module."""
        try:
            await self.bot.reload_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')

    _GIT_PULL_REGEX = re.compile(r'\s*(?P<filename>.+?)\s*\|\s*[0-9]+\s*[+-]+')

    def find_modules_from_git(self, output: str) -> list[tuple[int, str]]:
        files = self._GIT_PULL_REGEX.findall(output)
        ret: list[tuple[int, str]] = []
        for file in files:
            root, ext = os.path.splitext(file)
            if ext != '.py':
                continue

            if root.startswith('cogs/'):
                # A submodule is a directory inside the main cog directory for
                # my purposes
                ret.append((root.count('/') - 1, root.replace('/', '.')))

        # For reload order, the submodules should be reloaded first
        ret.sort(reverse=True)
        return ret

    async def reload_or_load_extension(self, module: str) -> None:
        try:
            await self.bot.reload_extension(module)
        except commands.ExtensionNotLoaded:
            await self.bot.load_extension(module)

    @_reload.command(name='all', hidden=True)
    async def _reload_all(self, ctx: commands.Context):
        """Reloads all modules, while pulling from git."""

        async with ctx.typing():
            stdout, stderr = await self.run_process('git pull')

        # progress and stuff is redirected to stderr in git pull
        # however, things like "fast forward" and files
        # along with the text "already up-to-date" are in stdout

        if stdout.startswith('Already up-to-date.'):
            return await ctx.send(stdout)

        modules = self.find_modules_from_git(stdout)
        mods_text = '\n'.join(f'{index}. `{module}`' for index, (_, module) in enumerate(modules, start=1))
        prompt_text = f'This will update the following modules, are you sure?\n{mods_text}'
        confirm = await ctx.prompt(prompt_text)
        if not confirm:
            return await ctx.send('Aborting.')

        statuses = []
        for is_submodule, module in modules:
            if is_submodule:
                try:
                    actual_module = sys.modules[module]
                except KeyError:
                    statuses.append((ctx.tick(None), module))
                else:
                    try:
                        importlib.reload(actual_module)
                    except Exception as e:
                        statuses.append((ctx.tick(False), module))
                    else:
                        statuses.append((ctx.tick(True), module))
            else:
                try:
                    await self.reload_or_load_extension(module)
                except commands.ExtensionError:
                    statuses.append((ctx.tick(False), module))
                else:
                    statuses.append((ctx.tick(True), module))

        await ctx.send('\n'.join(f'{status}: `{module}`' for status, module in statuses))

async def setup(bot):
    await bot.add_cog(Reload(bot))
