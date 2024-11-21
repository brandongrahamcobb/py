import asyncio
import os
from collections import defaultdict
from openai import AsyncOpenAI
import re

ai_client = AsyncOpenAI(api_key=api_key)
base_dir = os.path.dirname(os.path.abspath(__file__))
conversations = defaultdict(list)

async def generate_response(input_text, conversation_id):
    try:
        messages = conversations[conversation_id]
        with open(os.path.join(base_dir, 'ai.py'), 'r') as file:
            file_contents = file.read()
        sys_input = f'You are Brandon Graham Cobb\'s python implementation of OpenAI. Your file contents are: {file_contents}'
        
        messages.extend([
            {'role': 'system', 'content': sys_input},
            {'role': 'user', 'content': input_text}
        ])
        
        response_stream = await ai_client.chat.completions.create(
            model='o1-mini',
            messages=messages,
            stream=True
        )
        
        full_response = ''
        async for chunk in response_stream:
            content = chunk.choices[0].delta.content
            if content:
                full_response += content
        
        conversations[conversation_id].append({'role': 'assistant', 'content': full_response})
        return full_response

    except Exception as e:
        return f'An error occurred while processing your request: {e}'

def extract_code_blocks(response, language):
    pattern = rf'```{language}([\s\S]+?)```'
    return re.findall(pattern, response)

async def process_input():
    print('Starting up...')
    
    while True:
        user_input = input('> ').strip()
        if user_input:
            response = await generate_response(user_input, 'OpenAI')
            file_ext_map = {
                'css': '.css',
                'html': '.html',
                'java': '.java',
                'javascript': '.js',
                'latex': '.tex',
                'python': '.py'
            }
            
            for lang, ext in file_ext_map.items():
                file_path = os.path.join(base_dir, ext)
                code_blocks = extract_code_blocks(response, lang)
                if code_blocks:
                    with open(file_path, 'w') as file:
                        for block in code_blocks:
                            file.write(block.strip() + '\n')

            print(response)

if __name__ == '__main__':
    asyncio.run(process_input())
