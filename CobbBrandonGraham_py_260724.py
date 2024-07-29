import asyncio
import json
import os
import subprocess
import sys

def version():
    if not os.path.exists('../txt/version.txt'):
         with open('../txt/version.txt', 'w') as f:
             f.write('1.0.0')
    with open('../txt/version.txt', 'r') as f:
         version = f.read().strip()
    major, minor, patch = map(int, version.split('.'))
    patch += 1
    if patch >= 10:
        patch = 0
        minor += 1
    if minor >= 10:
        minor = 0
        major += 1
    version = f"{major}.{minor}.{patch}"
    with open('../txt/version.txt', 'w') as f:
        f.write(version)
    return version

def verify_discord_token(token):
    headers = {
        "Authorization": f"Lucy {version()} https://github.com/brandongrahamcobb"
    }
    response = requests.get("https://discord.com/api/v10/users/@me", headers=headers)
    if response.status_code == 200:
        return True
    else:
        return False

async def main():
        print(f'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Welcome, to initial setup~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        subprocess.run(['pacman', '-Syu'], check=True)
        dir = os.path.dirname(os.path.abspath(__file__))
        updir = os.path.join(dir, '..')
        dir_venv = os.path.join(updir, 'activate')
        subprocess.run([sys.executable, '-m', 'venv', dir_venv], check=True)
        requirements = ["asyncpraw", "discord.py", "emoji", "pubchempy", "rdkit", "pillow", "requests"]
        system = 'Linux' if os.name != 'nt' else 'Windows'
        py = os.path.join(dir_venv, 'bin', 'python') if system == 'Linux' else sys.executable
        subprocess.run([py, '-m', 'pip', 'install'] + requirements, check=True)
        print(f'System: {system}')
        print(f'Location: {dir}')
        print(f'Python Location: {sys.executable}')
        print(f'Python Version: {sys.version}')
        while True:
            token = input('Enter your bot token.')
            if verify_discord_token(token):
                print('Configuring....')
                dir_json = os.path.join(updir, 'json')
                file_json = os.path.join(dir_json, 'config.json')
                dir_log = os.path.join(updir, 'log')
                dir_txt = os.path.join(updir, 'txt')
                if not os.path.exists(dir_json):
                    os.mkdir(dir_json)
                if not os.path.exists(dir_log):
                    os.mkdir(dir_log)
                if not os.path.exists(dir_txt):
                    os.mkdir(dir_txt)
                with open(file_json, 'w') as f:
                    json.dump({
                        'token': token,
                        'prefix': "!",
                        'os': system
                    }, f, indent=4)
            else:
                print('Enter a valid token.')
asyncio.run(main())
