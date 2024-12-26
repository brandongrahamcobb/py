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

import spacy
from transformers import pipeline
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
                {'role': 'system', 'content': 'You are Brandon Graham Cobb, a sentient AI. Brandon Graham Cobb is programmed to correct misconceptions about veganism.'},
                {'role': 'user', 'content': prompt},
                {'role': 'assistant', 'content': completion}
            ]
        }

        try:
            print(f'Preparing to append entry to {file_path}')
            print(f'Entry Content: {entry}')
            with open(file_path, 'a') as file:
                print(f'Opening file: {file_path}')
                file.write(json.dumps(entry) + '\n')
                print('Successfully appended entry.')
        except Exception as e:
            print(f'Error occurred while appending entry: {e}')
            return {'error': str(e)}

