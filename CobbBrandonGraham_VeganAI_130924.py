from collections import defaultdict
from openai import AsyncOpenAI
from quart import Quart, request, Response, render_template
from uuid import uuid4

import asyncio
import os

ai_client = AsyncOpenAI(api_key='sk-proj-XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
app = Quart(__name__)
base = os.path.dirname(os.path.abspath(__file__))
conversations = defaultdict(list)

async def ai(input, conversation_id):
    try:
        messages = conversations[conversation_id]
        messages.append({'role': 'system', 'content': 'You are Brandon Graham Cobb\'s python implementation of OpenAI.'})
        messages.append({'role': 'user', 'content': input})
        stream = await ai_client.chat.completions.create(
            model='gpt-4o-mini',
            messages=messages,
            stream=True
        )
        full_response = ''
        async for chunk in stream:
            content = chunk.choices[0].delta.content
            if content is not None:
                full_response += content
        conversations[conversation_id].append({'role': 'assistant', 'content': full_response})
        parts = full_response.split('```')
        return parts
    except Exception as e:
        return 'An error occurred while processing your request.'

@app.route('/')
async def index():
    return await render_template('index.html')

@app.route('/api/main', methods=['POST'])
async def main():
    input = (await request.get_json())['input']
    conversation_id = uuid4()
    conversations[conversation_id].append({'role': 'user', 'content': input})
    async def event_stream():
        parts = await ai(input, conversation_id)
        for part in parts:
            yield part
    return Response(event_stream(), mimetype='text/event-stream')

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
