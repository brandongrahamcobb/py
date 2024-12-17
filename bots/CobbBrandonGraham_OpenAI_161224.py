#!/usr/bin/env python3
# master_script.py

import sys
import asyncio
import traceback
import uuid
import json
import os
from collections import defaultdict
from datetime import datetime as dt
from os.path import join
import yaml
import aiohttp

from PyQt5.QtWidgets import (
    QApplication,
    QWidget,
    QVBoxLayout,
    QTextEdit,
    QLineEdit,
    QPushButton,
    QLabel,
    QMessageBox,
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal

# ==========================
# Configuration and Helpers
# ==========================

def load_yaml(file_path):
    """
    Load YAML configuration from the given file path.
    """
    if os.path.exists(file_path):
        with open(file_path, 'r') as f:
            try:
                data = yaml.safe_load(f)
                return data
            except yaml.YAMLError as e:
                print(f"Error loading YAML file: {e}")
                return None
    else:
        print(f"Configuration file not found at {file_path}")
        return None

# ==========================
# API Request Handler Thread
# ==========================

class APIRequestThread(QThread):
    """
    QThread to handle asynchronous API requests without blocking the GUI.
    """
    response_received = pyqtSignal(str)
    error_occurred = pyqtSignal(str)

    def __init__(self, request_data, headers, url, save_path):
        super().__init__()
        self.request_data = request_data
        self.headers = headers
        self.url = url
        self.save_path = save_path

    async def make_request(self):
        try:
            async with aiohttp.ClientSession() as session:
                async with session.post(
                    url=self.url,
                    headers=self.headers,
                    json=self.request_data
                ) as response:
                    if response.status == 200:
                        result = await response.json()
                        full_response = ""
                        if "choices" in result:
                            for choice in result["choices"]:
                                content = choice.get("message", {}).get("content", "")
                                full_response += content
                        else:
                            self.error_occurred.emit("Invalid response structure from API.")
                            return
                        # Emit the response to the main thread
                        self.response_received.emit(full_response.strip())
                        # Save the response to a file
                        with open(self.save_path, 'a') as f:
                            f.write(full_response.strip() + "\n\n")
                    else:
                        error_text = await response.text()
                        self.error_occurred.emit(f"Error {response.status}: {error_text}")
        except Exception as e:
            tb = traceback.format_exc()
            self.error_occurred.emit(f"Exception: {str(e)}\n{tb}")

    def run(self):
        asyncio.run(self.make_request())

# ==========================
# Main GUI Class
# ==========================

class ChatWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("OpenAI Chat Interface")
        self.setGeometry(100, 100, 600, 500)

        self.layout = QVBoxLayout()

        # Label for input
        self.input_label = QLabel("Enter your message:")
        self.layout.addWidget(self.input_label)

        # Text input
        self.input_text = QLineEdit()
        self.layout.addWidget(self.input_text)

        # Send button
        self.send_button = QPushButton("Send")
        self.send_button.clicked.connect(self.handle_send)
        self.layout.addWidget(self.send_button)

        # Response display
        self.response_display = QTextEdit()
        self.response_display.setReadOnly(True)
        self.layout.addWidget(self.response_display)

        self.setLayout(self.layout)

        # Load configuration
        self.config = self.load_config()
        if not self.config:
            QMessageBox.critical(self, "Configuration Error", "Failed to load configuration.")
            sys.exit(1)

        # Set up API parameters
        self.setup_api_parameters()

        # Initialize conversation history
        self.conversations = defaultdict(list)
        self.custom_id = uuid.uuid4().hex

        # Response file path
        self.path_responses = join(self.config.get('path_home', os.path.expanduser('~')), 'response.txt')

    def load_config(self):
        """
        Load the YAML configuration file.
        """
        path_home = os.path.expanduser('~')
        path_config_yaml = os.path.join(path_home, '.config', 'openai', 'config.yaml')
        config = load_yaml(path_config_yaml)
        if not config:
            print("Failed to load configuration.")
        return config

    def setup_api_parameters(self):
        """
        Set up API request parameters based on the configuration.
        """
        api_keys = self.config.get('api_keys', {})
        self.API_KEY = api_keys.get('api_key_1')
        if not self.API_KEY:
            QMessageBox.critical(self, "API Key Error", "API key not found in configuration.")
            sys.exit(1)

        self.headers = {
            'Authorization': f'Bearer {self.API_KEY}',
            'Content-Type': 'application/json',
            'OpenAI-Project': self.config.get('openai_project', 'default_project'),
            'User-Agent': self.config.get('user_agent', 'user@example.com')
        }

        # Model settings
        self.MODEL = self.config.get('model', 'gpt-4o-mini')

        self.MODEL_OUTPUT_LIMITS = {
            "gpt-3.5-turbo": 4096,
            "gpt-4": 8192,
            "gpt-4-32k": 32768,
            "gpt-4o": 4096,
            "gpt-4o-mini": 16384,
            "gpt-4-turbo": 4096,
            "o1-preview": 32768,
            "o1-mini": 65536,
        }

        self.MODELS_CURRENT = ['o1-preview', 'o1-mini']
        self.MODELS_DEPRECATED = [
            'gpt-3.5-turbo', 'gpt-4', 'gpt-4-32k',
            'gpt-4o', 'gpt-4o-mini', 'gpt-4-turbo'
        ]

        self.MAX_TOKENS = self.MODEL_OUTPUT_LIMITS.get(self.MODEL, 4096)
        self.URL = self.config.get('api_url', "https://api.openai.com/v1/chat/completions")

    def handle_send(self):
        """
        Handle the send button click event.
        """
        input_text = self.input_text.text().strip()
        if not input_text:
            QMessageBox.warning(self, "Input Error", "Please enter a message.")
            return

        # Display user's message
        self.display_message("User", input_text)

        # Clear input field
        self.input_text.clear()

        # Prepare request data
        request_data = {
            "metadata": {
                "user": self.config.get('user_name', 'User'),
                "timestamp": str(dt.utcnow())
            },
            "model": self.MODEL,
            "n": 1,
            "stop": None,
            "store": True,
            "stream": False,
            "temperature": 0.7,
            "top_p": 1.0,
        }

        if self.MODEL in self.MODELS_CURRENT:
            request_data.update({
                "max_completion_tokens": self.MAX_TOKENS,
                'messages': [{'role': 'user', 'content': input_text}]
            })
        else:
            request_data.update({"max_tokens": self.MAX_TOKENS})

        # If system input is defined
        sys_input = self.config.get('sys_input')
        if sys_input:
            request_data.update({
                'messages': [
                    {'role': 'user', 'content': input_text},
                    {'role': 'system', 'content': sys_input}
                ]
            })
        else:
            request_data.update({
                'messages': [{'role': 'user', 'content': input_text}]
            })

        # Start the API request thread
        self.api_thread = APIRequestThread(
            request_data=request_data,
            headers=self.headers,
            url=self.URL,
            save_path=self.path_responses
        )
        self.api_thread.response_received.connect(self.handle_response)
        self.api_thread.error_occurred.connect(self.handle_error)
        self.api_thread.start()

    def handle_response(self, response):
        """
        Handle the API response.
        """
        self.display_message("Assistant", response)

    def handle_error(self, error_message):
        """
        Handle errors during the API request.
        """
        QMessageBox.critical(self, "API Error", error_message)

    def display_message(self, role, message):
        """
        Display a message in the response display area.
        """
        if role == "User":
            self.conversations[self.custom_id].append({'role': 'user', 'content': message})
            display_text = f"<b>User:</b> {message}\n"
        elif role == "Assistant":
            self.conversations[self.custom_id].append({'role': 'assistant', 'content': message})
            display_text = f"<b>Assistant:</b> {message}\n"
        else:
            display_text = f"<b>{role}:</b> {message}\n"

        self.response_display.append(display_text)
        self.response_display.moveCursor(Qt.TextCursor.End)

# ==========================
# Entry Point
# ==========================

def main():
    app = QApplication(sys.argv)
    window = ChatWindow()
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
