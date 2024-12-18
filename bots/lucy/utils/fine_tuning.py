import aiofiles
from datetime import datetime as dt
import json
import os
import traceback

async def generate_training_file(conversations):
    training_data = []
    placeholder = "carnist nonsense"
    for user_id, messages in conversations.items():
        previous_user_message = None  # Track the last user message.
        for message in messages:
            try:
                if isinstance(message, dict) and 'role' in message and 'content' in message:
                    if message['role'] == 'user':
                        previous_user_message = message['content']
                    elif message['role'] == 'assistant' and previous_user_message:
                        assistant_response = message['content']
                        training_entry = {
                            "prompt": placeholder,
                            "completion": f" {assistant_response}"
                        }
                        training_data.append(training_entry)
                        previous_user_message = None  # Reset after pairing
                else:
                    print(f"Invalid message format: {message}")
            except Exception as e:
                print(f"Error processing message {message}: {e}")
    try:
        home = os.path.expanduser('~')
        async with aiofiles.open(os.path.join(home, 'Downloads', f'training{dt.now()}.jsonl'), 'w', encoding='utf-8') as f:
            for entry in training_data:
                await f.write(json.dumps(entry) + '\n')
    except Exception as e:
        traceback.print_exc()
