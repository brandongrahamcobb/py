""" helpers.py
    Copyright (C) 2024 github.com/brandongrahamcobb

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
"""

import logging
import logging.handlers

def setup_logging():
    global logger
    log_file = 'discord.log'
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.addHandler(file_handler)
    if not os.path.exists(log_file):
        open(log_file, 'a').close()

import json
import os

def load_config():
    config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config.json')
    if not os.path.exists(config_path):
        raise FileNotFoundError("Configuration file not found.")
    with open(config_path, 'r') as f:
        return json.load(f)

import requests
from bs4 import BeautifulSoup

def fetch_and_parse(url):
    try:
        response = requests.get(url)
        response.raise_for_status()
        soup = BeautifulSoup(response.content, 'html.parser')
        return soup.get_text()
    except requests.RequestException as e:
        print(f"Request failed: {e}")
        return None

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
import time

def interact_with_chatgpt(prompt):
    driver = webdriver.Chrome()  # Or use another driver
    driver.get('https://chat.openai.com/')  # OpenAI ChatGPT URL
    time.sleep(10)  # Wait for page load
    input_field = driver.find_element(By.CSS_SELECTOR, 'textarea')
    input_field.send_keys(prompt + Keys.RETURN)
    time.sleep(5)  # Wait for response
    messages = driver.find_elements(By.CSS_SELECTOR, '.message')
    response = messages[-1].text if messages else 'No response found'
    driver.quit()
    return response
