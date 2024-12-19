from bs4 import BeautifulSoup
import utils.helpers as helpers

import json
import requests

def script(version: str, reference: str):
    BIBLE_IDS = {
         'esv': 'de4e12af7f28f599-02',
         'nkjv': 'de4e12af7f28f599-01',
         'niv': '06125adad2d5898a-01',
    }
    version = version.lower()
    if version in BIBLE_IDS:
        bible_id = BIBLE_IDS[version]
        api = f'https://api.scripture.api.bible/v1/bibles/{bible_id}/search?query={reference}'
        response = requests.get(api, headers=helpers.SCRIPTURE_HEADERS)
        if response.ok:
            json = response.json()
            passages = json.get('data', {}).get('passages', [])
            soup = BeautifulSoup(passages[0].get('content'), 'html.parser')
            soup.get_text()
            cleaned_content = soup.get_text()
            message = f'**{reference}** ({version.upper()})\n{cleaned_content}'
            return message
        if version == '':
            response = requests.get(f'https://api.alquran.cloud/v1/ayah/{reference}/en.asad', headers=get_scripture_headers())
            if response.ok:
                json = response.json()
                message = f'**{reference}** ({version.upper()})\n{json['data']['text']}'
                return message

