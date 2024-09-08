from openai import AsyncOpenAI
import openai
import os
from quart import Quart, request, Response, render_template
import asyncio
import requests
import logging

app = Quart(__name__)
openai_client = AsyncOpenAI()

@app.route('/')
async def index():
    return await render_template('index.html')

logging.basicConfig(level=logging.INFO)

@app.route('/api/gpt', methods=['POST'])
async def chatgpt():
    logging.info("Received request for GPT response.")
    user_input = (await request.get_json())['input']
    logging.info(f"User input: {user_input}")
    async def generate_response():
        try:
            stream = await openai_client.chat.completions.create(
                model='gpt-4o-mini',
                messages=[
                    {'role': 'system', 'content': 'You are vegan.'},
                    {'role': 'user', 'content': user_input}
                ],
                stream=True
            )
            async for chunk in stream:
                if chunk.choices[0].delta.content is not None:
                    print(chunk.choices[0].delta.content)
                    yield chunk.choices[0].delta.content
        except Exception as e:
            logging.error(f"Error occurred: {e}")
            yield "An error occurred while processing your request."
    return Response(generate_response(), mimetype='text/event-stream')

if __name__ == '__main__':
    app.run(debug=True)
