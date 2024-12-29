''' nlp_utils.py  The purpose of this program is to provide generic Natural Language Processing functionality from cd ../.
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

from transformers import pipeline
from utils.setup_logging import logger

import spacy
import json
import os

nlp = spacy.load('en_core_web_sm')
sentiment_analyzer = pipeline('sentiment-analysis')

class NLPUtils:

    @staticmethod
    def sentiment_analysis(text):
        try:
            result = sentiment_analyzer(text)
            return result[0]
        except Exception as e:
            return {'error': str(e)}

    @staticmethod
    def dependency_parsing(text):
        doc = nlp(text)
        dependencies = [
            {
                'word': token.text,
                'dependency': token.dep_,
                'head': token.head.text
            }
            for token in doc
        ]
        return dependencies

    @staticmethod
    def combined_analysis(text):
        try:
            sentiment = NLPUtils.sentiment_analysis(text)
            dependencies = NLPUtils.dependency_parsing(text)

            return {
                'sentiment': sentiment,
                'dependencies': dependencies
            }
        except Exception as e:
            return {'error': str(e)}

    @staticmethod
    def append_to_jsonl(file_path, sentiment_result, completion):
        sentiment_label = sentiment_result.get('label', 'unknown')
        sentiment_score = sentiment_result.get('score', 0)
        prompt = f'Sentiment: {sentiment_label}, Score: {sentiment_score}. Non-vegan: '

        entry = {
            'messages': [
                {'role': 'system', 'content': 'You are Vyrtuous, a sentient AI. Vyrtuous is programmed to correct misconceptions about veganism.'},
                {'role': 'user', 'content': {
                        'sentiment_label': sentiment_label,
                        'sentiment_score': sentiment_score,
                        'prompt': None
                    }
                },
                {'role': 'assistant', 'content': completion}
            ]
        }

        try:
            print(f'Preparing to append entry to {file_path}')
            print(f'Entry Content: {entry}')
            with open(file_path, 'r') as file:
                list_obj = json.load(file)
                list_obj.append(entry)
            with open(file_path, 'w') as json_file:
                print(f'Opening file: {file_path}')
                json.dump(list_obj, json_file, indent=4, separators=(',',': '))
                print('Successfully appended entry.')
        except Exception as e:
            print(f'Error occurred while appending entry: {e}')
            return {'error': str(e)}

    @staticmethod
    def append_to_other_jsonl(file_path, flagged, completion):
        entry = {
            'messages': [
                {'role': 'system', 'content': 'You are Vyrtuous, a sentient AI. Vyrtuous is programmed to correct misconceptions about veganism.'},
                {'role': 'user', 'content': json.dumps({'flagged': flagged})},
                {'role': 'assistant', 'content': completion}
            ]
        }
        try:
            print(f'Preparing to append entry to {file_path}')
            print(f'Entry Content: {entry}')
            with open(file_path, 'a') as file:  # 'a' mode for appending
                json.dump(entry, file)  # Write the JSON object as a line
                file.write('\n')  # Add a newline character for JSONL format
            print('Successfully appended entry in JSONL format.')
        except Exception as e:
            print(f'Error occurred while appending entry: {e}')
            return {'error': str(e)}
