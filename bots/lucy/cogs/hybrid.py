''' hybrid.py The purpose of this program is to provide the core functionality to Vyrtuous.
    Copyright (C) 2024  github.com/brandongrahamcobb

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
from discord.ext import commands, tasks
from discord.ext.commands import Parameter
from datetime import datetime
from dataclasses import dataclass, field
from os.path import abspath, dirname, expanduser, join

import uuid
import aiohttp
import asyncio
import bot.utils.api_request_parallel_processor as api_processor
import bot.utils.openai_helpers as openai_helpers
import bot.utils.helpers as helpers
import discord
import json
import os
import time
import logging
import re

CHECKMARK_EMOJI = '✅'
CROSS_EMOJI = '❌'
DEFAULT_PREFS = {
    "completions": 1,
    "conversation_id": "",
    "max_tokens": 150,
    "model": "gpt-3.5-turbo",
    "stop": None,
    "store": False,
    "stream": True,
    "sys_input": None,
    "temperature": 0.7,
    "top_p": 1.0,
}
dir_base = dirname(abspath(__file__))
path_helpers = join(dir_base, '..', 'utils', 'openai_helpers.py')
path_home = expanduser('~')
path_openai_helpers = join(dir_base, '..', 'utils', 'openai_helpers.py')
path_hybrid = join(dir_base, '..', 'cogs', 'hybrid.py')
path_indica = join(dir_base, '..', 'cogs', 'indica.py')
path_main = join(dir_base, '..', 'main.py')
path_sativa = join(dir_base, '..', 'cogs', 'sativa.py')

class Hybrid(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.status_tracker=self.bot.status_tracker
        self.config = self.bot.config
        self.indica = self.bot.get_cog('Indica')
        self.path_requests = join(dir_base, '..',  'requests.jsonl')
        self.path_results = join(dir_base, '..',  'results.jsonl')

        self.sativa = self.bot.get_cog('Sativa')
        self.helpers_py = helpers.load_contents(path_helpers)
        self.openai_helpers_py = helpers.load_contents(path_openai_helpers)
        self.hybrid_py = helpers.load_contents(path_hybrid)
        self.indica_py = helpers.load_contents(path_indica)
        self.main_py = helpers.load_contents(path_main)
        self.sativa_py = helpers.load_contents(path_sativa)
        self.sys_input = f"""
            Your main.py file is {self.main_py}.
            Your cogs are in cogs/ {self.hybrid_py}, {self.indica_py}, {self.sativa_py}.
            Your helpers are in utils/ {self.helpers_py}.
            Your openai_helpers are in utils/ {self.openai_helpers_py}.
        """
        self.requests_file = join(dir_base, '..',  'requests.jsonl')
        self.results_file = join(dir_base, '..',  'results.jsonl')
        # Replace with your API headers and endpoint
        self.api_key = self.config['api_keys']['api_key_1']
        self.api_url = 'https://api.openai.com/v1/chat/completions'
        self.status_tracker = api_processor.StatusTracker()
        self.processing_task = None

    @commands.hybrid_command(name="batch", help="Processes a batch from requests.jsonl file.")
    async def batch(self, ctx: commands.Context):
        """Processes a batch using OpenAI's Batch API."""
        await ctx.defer()  # Defer response for hybrid commands
        success, message = await openai_helpers.process_batch_from_file(self.path_requests)

        if success:
            await ctx.send(f"Batch processed successfully:\n{message}")
        else:
            await ctx.send(f"Batch processing failed: {message}")

    @commands.hybrid_command(name="process")
    async def start_processing(self, ctx: commands.Context):
        if "processing_task" in self.bot.__dict__ and not self.bot.processing_task.done():
            await ctx.send("Processing is already running.")
            return
        self.bot.processing_task = asyncio.create_task(
            api_processor.process_api_requests_from_file(
                requests_filepath=self.requests_file,
                save_filepath=self.results_file,
                request_url=self.api_url,
                api_key=self.api_key,
                max_requests_per_minute=60,
                max_tokens_per_minute=100000,
                max_attempts=3,
                token_encoding_name="cl100k_base",
                logging_level=1,
                status_tracker=self.status_tracker,
            )
        )
        await ctx.send("Processing has started!")

    @commands.hybrid_command(name='request')
    async def request_openai(
        self,
        ctx: commands.Context,
        input_text: str = commands.parameter(default='', description='Input text'),
        completions: int = commands.parameter(default=1, description='Number of completions'),
        conversation_id: str = commands.parameter(default='', description='Conversation id'),
        max_tokens: int = commands.parameter(default=150, description='Max tokens'),
        model: str = commands.parameter(default='gpt-3.5-turbo', description="""
            gpt-4o-mini, gpt-4, gpt-4o
        """),
        stop: str = commands.parameter(default=None, description='Stop'),
        store: bool = commands.parameter(default=False, description='Queue or request completions'),
        stream: bool = commands.parameter(default=False, description='Stream On/Off'),
        sys_input: str = commands.parameter(default=None, description='System input'),
        temperature: float = commands.parameter(default=0.7, description='Temperature'),
        top_p: float = commands.parameter(default=1.0, description='Nucleus Sampling')
    ):
        """Handles a request to the OpenAI GPT model."""
    
        # Convert Parameter objects to their defaults if necessary.
        # (Only needed if you're still having Parameter objects issue)
        from discord.ext.commands import Parameter
        def param_to_value(p): return p.default if isinstance(p, Parameter) else p
        input_text = param_to_value(input_text)
        completions = param_to_value(completions)
        conversation_id = param_to_value(conversation_id)
        max_tokens = param_to_value(max_tokens)
        model = param_to_value(model)
        stop = param_to_value(stop)
        store = param_to_value(store)
        stream = param_to_value(stream)
        sys_input = param_to_value(sys_input)
        temperature = param_to_value(temperature)
        top_p = param_to_value(top_p)
    
        # Load user preferences
        prefs = helpers.load_json(helpers.path_preferences).get(str(ctx.author.id), DEFAULT_PREFS)
    
        # Resolve final parameters
        final_completions = prefs['completions'] if prefs.get('completions') is not None else completions
        final_conversation_id = prefs['conversation_id'] if prefs.get('conversation_id') is not None else conversation_id
        final_max_tokens = prefs['max_tokens'] if prefs.get('max_tokens') is not None else max_tokens
        final_model = prefs['model'] if prefs.get('model') is not None else model
        final_stop = prefs['stop'] if prefs.get('stop') is not None else stop
        final_store = prefs['store'] if prefs.get('store') is not None else store
        final_stream = prefs['stream'] if prefs.get('stream') is not None else stream
        final_sys_input = prefs['sys_input'] if prefs.get('sys_input') is not None else sys_input
        final_temperature = prefs['temperature'] if prefs.get('temperature') is not None else temperature
        final_top_p = prefs['top_p'] if prefs.get('top_p') is not None else top_p

        # Define which models can have system input
        MODELS_WITH_SYS_INPUT = {"gpt-3.5-turbo", "gpt-4-turbo", "gpt-4o", "gpt-4o-mini"}
    
        # Construct the messages list
        messages = [{"role": "user", "content": input_text}]
    
        # Only include a system message if the model supports it
        if final_model in MODELS_WITH_SYS_INPUT:
            if final_sys_input is None:
                final_sys_input = "You are a helpful assistant."
            # Insert the system message at the beginning
            messages.insert(0, {"role": "system", "content": final_sys_input})
    
        request_data = {
            "max_tokens": final_max_tokens,
            "messages": messages,
            "model": final_model,
            "temperature": final_temperature,
            "top_p": final_top_p,
            "n": final_completions,
            "stop": final_stop,
            "store": final_store,
            "stream": final_stream,
        }
    
        # Store the request if batching is required
        if final_store:
            request_data.update({
                "custom_id": f"{ctx.author.id}-{uuid.uuid4().hex}",
                "method": "POST",
                "url": "/v1/chat/completions",
                "metadata": {"user": str(ctx.author.id), "timestamp": str(datetime.utcnow())}
            })
            api_processor.append_to_jsonl(request_data, self.path_requests)
            await ctx.send(f"Your request has been queued for batch processing with ID `{request_data['custom_id']}`.")
            return
    
        # Immediate processing
        async for response in openai_helpers.create_completion(
            completions=final_completions,
            conversation_id=final_conversation_id,
            input_text=input_text,
            max_tokens=final_max_tokens,
            model=final_model,
            stop=final_stop,
            store=final_store,
            stream=final_stream,
            sys_input=final_sys_input if final_model in {"gpt-3.5-turbo", "gpt-4", "gpt-4-turbo"} else '',
            temperature=final_temperature,
            top_p=final_top_p,
        ):
            message = await ctx.send(response.strip())
            await message.add_reaction("✅")
            await message.add_reaction("❌")
    
    @commands.hybrid_command(name='set')
    async def set_openai(
        self,
        ctx: commands.Context,
        completions: int = commands.parameter(default=1, description='Number of completions'),
        conversation_id: str = commands.parameter(default='', description='Conversation ID'),
        max_tokens: int = commands.parameter(default=150, description='Max tokens'),
        model: str = commands.parameter(default='gpt-3.5-turbo', description='Model'),
        stop: str = commands.parameter(default=None, description='Stopping criteria.'),
        store: bool = commands.parameter(default=False, description='Queue or request completions'),
        stream: bool = commands.parameter(default=True, description='Stream On/Off'),
        sys_input: str = commands.parameter(default=None, description='System input'),
        temperature: float = commands.parameter(default=0.7, description='Temperature'),
        top_p: float = commands.parameter(default=1.0, description='Nucleus sampling')
    ):
        preferences = {
            "completions": completions,
            "conversation_id": conversation_id,
            "max_tokens": max_tokens,
            "model": model,
            "stop": stop,
            "store": store,
            "stream": stream,
            "sys_input": sys_input,
            "temperature": temperature,
            "top_p": top_p,
        }
        prefs = helpers.load_json(helpers.path_preferences) or {}
        prefs[str(ctx.author.id)] = preferences
        await helpers.save_json(helpers.path_preferences, prefs)
        await ctx.send("Your preferences have been updated!")

    @commands.hybrid_command(name="status")
    async def status(self, ctx):
        status_msg = (
            f"Tasks Started: {self.status_tracker.num_tasks_started}\n"
            f"Tasks In Progress: {self.status_tracker.num_tasks_in_progress}\n"
            f"Tasks Succeeded: {self.status_tracker.num_tasks_succeeded}\n"
            f"Tasks Failed: {self.status_tracker.num_tasks_failed}\n"
            f"API Errors: {self.status_tracker.num_api_errors}\n"
            f"Rate Limit Errors: {self.status_tracker.num_rate_limit_errors}\n"
            f"Other Errors: {self.status_tracker.num_other_errors}"
        )
        await ctx.send(f"**Processing Status:**\n{status_msg}")

async def setup(bot: commands.Bot):
    await bot.add_cog(Hybrid(bot))
