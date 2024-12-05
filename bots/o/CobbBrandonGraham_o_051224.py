import sys
import asyncio
import traceback
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QTextEdit, QLineEdit
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from openai import AsyncOpenAI
import openai
import os
import yaml

# Global conversations dictionary
conversations = {}

path_home = os.path.expanduser('~')
path_config_yaml = os.path.join(path_home, '.config', 'vyrtuous', 'config.yaml')

def load_config():
    if os.path.exists(path_config_yaml):
        with open(path_config_yaml, 'r') as f:
            data = yaml.safe_load(f)
            if data:
                return data
    else:
        return None

async def create_completion(input_text, conversation_id):
    try:
        config = load_config()
        api_key = config['api_keys']['api_key_2']
        ai_client = AsyncOpenAI(api_key=api_key)
        messages = conversations.get(conversation_id, [])
        messages.append({'role': 'user', 'content': input_text})
        # Use the appropriate model here ('o1-mini' as per your requirement)
        stream = await ai_client.chat.completions.create(
            model='o1-mini',
            messages=messages,
            stream=True
        )
        full_response = ''
        async for chunk in stream:
            content = chunk.choices[0].delta.content
            if content is not None:
                full_response += content
                yield content
        messages.append({'role': 'assistant', 'content': full_response})
        conversations[conversation_id] = messages
    except Exception as e:
        yield traceback.format_exc()

async def deprecated_create_completion(input_text, sys_input, conversation_id):
    try:
        config = load_config()
        api_key = config['api_keys']['api_key_2']
        ai_client = AsyncOpenAI(api_key=api_key)
        messages = conversations.get(conversation_id, [])
        messages.append({'role': 'system', 'content': sys_input})
        messages.append({'role': 'user', 'content': input_text})
        # Use the appropriate model here ('gpt-4o-mini' as per your requirement)
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
                yield content
        messages.append({'role': 'assistant', 'content': full_response})
        conversations[conversation_id] = messages
    except Exception as e:
        yield traceback.format_exc()

class Worker(QThread):
    partial_result_signal = pyqtSignal(str)
    finished_signal = pyqtSignal()

    def __init__(self, input_text, conversation_id, sys_input=None, use_deprecated=False):
        super().__init__()
        self.input_text = input_text
        self.conversation_id = conversation_id
        self.sys_input = sys_input
        self.use_deprecated = use_deprecated

    def run(self):
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        if self.use_deprecated:
            loop.run_until_complete(self.run_deprecated_async())
        else:
            loop.run_until_complete(self.run_async())
        loop.close()

    async def run_async(self):
        try:
            async for content in create_completion(self.input_text, self.conversation_id):
                self.partial_result_signal.emit(content)
            self.finished_signal.emit()
        except Exception as e:
            self.partial_result_signal.emit(f"\nError: {str(e)}")
            self.finished_signal.emit()

    async def run_deprecated_async(self):
        try:
            async for content in deprecated_create_completion(self.input_text, self.sys_input, self.conversation_id):
                self.partial_result_signal.emit(content)
            self.finished_signal.emit()
        except Exception as e:
            self.partial_result_signal.emit(f"\nError: {str(e)}")
            self.finished_signal.emit()

class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Arch Linux System Manager with AI Assistant")
        self.layout = QVBoxLayout()

        # Text display area
        self.text_display = QTextEdit()
        self.text_display.setReadOnly(True)
        self.text_display.setAcceptRichText(False)

        # Input field
        self.input_field = QLineEdit()
        self.input_field.returnPressed.connect(self.handle_input)

        self.layout.addWidget(self.text_display)
        self.layout.addWidget(self.input_field)
        self.setLayout(self.layout)

    def handle_input(self):
        input_text = self.input_field.text()
        self.input_field.clear()

        # Decide which function to use based on input or other conditions
        # For example, if input_text starts with "use deprecated", use the deprecated function
        if input_text.startswith("@spawd"):
            self.use_deprecated = False
            sys_input = None
        else:
            self.use_deprecated = True
            sys_input = "You are Lucy, the vegan AI."
            # Remove the leading command
            input_text = input_text[len("!"):].strip()

        self.text_display.append(f'You: {input_text}')

        # Start the worker thread
        self.worker = Worker(
            input_text,
            'conversation_id_1',  # Use a unique conversation ID if needed
            sys_input=sys_input,
            use_deprecated=self.use_deprecated
        )
        self.worker.partial_result_signal.connect(self.display_partial_response)
        self.worker.finished_signal.connect(self.worker_finished)
        self.worker.start()

        # Prepare to display assistant's response
        self.text_display.append('Assistant: ')
        self.assistant_cursor = self.text_display.textCursor()
        self.assistant_cursor.movePosition(self.assistant_cursor.End)
        self.text_display.setTextCursor(self.assistant_cursor)

    def display_partial_response(self, content):
        self.assistant_cursor.insertText(content)
        self.text_display.setTextCursor(self.assistant_cursor)

    def worker_finished(self):
        # Optionally handle any cleanup after the worker finishes
        pass

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
