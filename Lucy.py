import subprocess
import venv
import json
import sys
import os
import requests

location = os.path.dirname(os.path.abspath(__file__))
requirements = ["discord.py", "pubchempy", "rdkit", "pillow", "requests"]
script_path = os.path.join(location, 'Main.py')
super_path = os.path.join(location, '..')
system = 'Linux' if os.name != 'nt' else 'Windows'
venv_path = os.path.join(location, 'venv')
venv_python = os.path.join(venv_path, 'bin', 'python') if system == 'Linux' else None

print(f'System: {system}')
print(f'Location: {location}')
print(f'Python Location: {sys.executable}')
print(f'Python Version: {sys.version}')


def clear():
    if system == 'Windows': os.system('cls')
    elif system == 'Linux': os.system('clear')

def verify_discord_token(token):
    headers = {
       'Authorization': f'Bot {token}'
    }
    url = 'https://discord.com/api/v9/users/@me'
    try:
        response = requests.get(url, headers=headers)
        response.raise_for_status()  # Raise an exception for 4xx or 5xx errors
        data = response.json()
        return data
    except requests.exceptions.HTTPError as err:
        print(f"HTTP error occurred: {err}")
        return None
    except Exception as err:
        print(f"Error occurred: {err}")
        return None

q = input('Is this initial setup? (y/n)')
if q == 'y':
    try:
        if system == 'Linux':
            if not os.path.exists(venv_path):
                q = input('Do you want to create a venv? (y/n): ')
                if q.lower() == 'y':
                    print('Creating virtual environment...')
                    venv.create(venv_path, with_pip=True)
                    print('Virtual environment created successfully.')
            print('Activating virtual environment...')
            subprocess.run([venv_python, '-m', 'venv', 'activate'], check=True)
            print('Installing dependencies...')
            subprocess.run([venv_python, '-m', 'pip', 'install'] + requirements, check=True)
        else:
            print('Installing dependencies...')
            subprocess.run([sys.executable, '-m', 'pip', 'install'] + requirements, check=True)
        print('Creating json folder...')
        os.mkdir(os.path.join(super_path, 'json'))
        print('Creating log folder...')
        os.mkdir(os.path.join(super_path, 'log'))
        print('Creating txt folder...')
        os.mkdir(os.path.join(super_path, 'txt'))
        print('Creating config.json...')
        while True:
            token = input('Enter your bot token: ')
            if verify_discord_token(token):
                break
            print('Please enter a valid token.')
        print('Saving config...')
        with open('../json/config.json', 'w') as f:
            json.dump({
                'token': token,
                'prefix': "!",
                'os': system
            }, f, indent=4)
        f.close()
        print('Done configuring...')
    except Exception as e:
        None

if __name__ == '__main__':
    try:
        if system == 'Linux':
            subprocess.run([venv_python, script_path], check=True)
        else:
            subprocess.run([sys.executable, script_path], check=True)
    except Exception as e:
        None
