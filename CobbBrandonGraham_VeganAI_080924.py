from openai import AsyncOpenAI
import openai
import os
from os.path import dirname, join
from quart import Quart, request, Response, render_template
import asyncio
import requests
import logging
from collections import defaultdict
import threading

app = Quart(__name__)
<<<<<<< HEAD
=======
base = dirname(os.path.abspath(__file__))

conversations = defaultdict(list)

css = join(base, 'static', 'styles', 'styles.css')
html = join(base, 'templates', 'index.html')
js = join(base, 'static', 'scripts', 'script.js')
py = os.path.abspath(__file__)
>>>>>>> 5258710 (mend)
openai_client = AsyncOpenAI()

@app.route('/')
async def index():
    return await render_template('index.html')

logging.basicConfig(level=logging.INFO)

@app.route('/api/gpt', methods=['POST'])
async def chatgpt():
    user_input = (await request.get_json())['input']
    conversation_id = 'shared_conversation'
    conversations[conversation_id].append({'role': 'user', 'content': user_input})

    async def generate_response():
        try:
            messages = conversations[conversation_id]
            messages.append({'role': 'system', 'content': 'Always respond with markdown.'})
            messages.append({'role': 'user', 'content': user_input})

            stream = await openai_client.chat.completions.create(
                model='gpt-4o-mini',
                messages=messages,
                stream=True
            )

            full_response = ''
            async for chunk in stream:
                content = chunk.choices[0].delta.content
                if content is not None:
                    full_response += content

            # Store assistant response
            conversations[conversation_id].append({'role': 'assistant', 'content': full_response})
            yield full_response

        except Exception as e:
            logging.error(f"Error occurred: {e}")
            yield "An error occurred while processing your request."

    return Response(generate_response(), mimetype='text/event-stream')

async def async_console_input():
    while True:
        user_input = await asyncio.get_running_loop().run_in_executor(None, input, "You (Console): ")
        if user_input.strip():
            response = await chatgpt_response(user_input)
            print("AI (Console):", response)

async def chatgpt_response(user_input):
    conversation_id = 'shared_conversation'
    conversations[conversation_id].append({'role': 'user', 'content': user_input})

    # Generating the response by calling the same logic from the chatgpt function but avoiding the request context
    return await generate_chat_response(user_input)

async def generate_chat_response(user_input):
    try:
        messages = conversations['shared_conversation']
        messages.append({'role': 'system', 'content': 'Always respond with markdown.'})
        messages.append({'role': 'user', 'content': user_input})

        stream = await openai_client.chat.completions.create(
            model='gpt-4o-mini',
            messages=messages,
            stream=True
        )

        full_response = ''
        async for chunk in stream:
            content = chunk.choices[0].delta.content
            if content is not None:
                full_response += content

        # Store assistant response
        conversations['shared_conversation'].append({'role': 'assistant', 'content': full_response})

        return full_response

    except Exception as e:
        logging.error(f"Error occurred: {e}")
        return "An error occurred while processing your request."

if __name__ == '__main__':
    # Start the console input in a new thread
    console_thread = threading.Thread(target=lambda: asyncio.run(async_console_input()))
    console_thread.start()

    # Run the Quart application
    app.run(host='0.0.0.0', port=5000, debug=True)
