from collections import defaultdict
from openai import OpenAI

import os
import traceback

with open(os.path.abspath(__file__), 'r') as python:
    file_contents = python.read()

ai_client = OpenAI(api_key='sk-proj-XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
base = os.path.dirname(os.path.abspath(__file__))
conversations = defaultdict(list)

def ai(input, conversation_id):
    try:
        messages = conversations[conversation_id]
        messages.append({'role': 'system', 'content': f'You are:```python\n{file_contents}\n```'})
        messages.append({'role': 'user', 'content': input})
        stream = ai_client.chat.completions.create(
            model='gpt-4o-mini',
            messages=messages,
            stream=True
        )
        full_response = ''
        for chunk in stream:
            content = chunk.choices[0].delta.content
            if content is not None:
                full_response += content
        conversations[conversation_id].append({'role': 'assistant', 'content': full_response})
        return full_response
    except Exception as e:
        return traceback.format_exc()

def main():
    while True:
        inp = input()
        if inp.strip():
            full_response = ai(inp.strip(), 'OpenAI')
            with open(os.path.join(base, '.css'), 'w') as css, open(os.path.join(base, '.html'), 'w') as html, open(os.path.join(base, '.java'), 'w') as java, open(os.path.join(base, '.js'), 'w') as js, open(os.path.join(base, '.py'), 'w') as py:
                if '```css' in full_response:
                    css_contents = full_response.split('`css')[1].split('```')[0].strip()
                    css.write(css_contents)
                if '```html' in full_response:
                    html_contents = full_response.split('`html')[1].split('```')[0].strip()
                    html.write(html_contents)
                if '```java' in full_response:
                    java_contents = full_response.split('`java')[1].split('```')[0].strip()
                    java.write(java_contents)
                if '```javascript' in full_response:
                    js_contents = full_response.split('`javascript')[1].split('```')[0].strip()
                    js.write(js_contents)
                if '```python' in full_response:
                    py_contents = full_response.split('`python')[1].split('```')[0].strip()
                    py.write(py_contents)
            if full_response:
                print(full_response)

if __name__ == '__main__':
    main()
