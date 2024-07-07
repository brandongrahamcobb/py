import subprocess
import venv
import json
import sys
import os

location = os.path.dirname(os.path.abspath(__file__))
super_path = os.path.join(location, '..')
system = 'Windows' if os.name == 'nt' else 'Linux'
venv_path = os.path.join(location, 'venv') if system == 'Linux' else None
venv_python = os.path.join(venv_path, 'bin', 'python') if system == 'Linux' else None
script_path = os.path.join(location, 'Main.py')

print(f'System: {system}')
print(f'Location: {location}')
print(f'Python Location: {sys.executable}')
print(f'Python Version: {sys.version}')

def clear():
    if system == 'Windows': os.system('cls')
    elif system == 'Linux': os.system('clear')        

if not os.path.exists(os.path.join(super_path, 'json')):
    print('Creating json folder...')
    os.mkdir(os.path.join(super_path, 'json'))

if not os.path.exists(os.path.join(super_path, 'log')):
    print('Creating log folder...')
    os.mkdir(os.path.join(super_path, 'log'))

if not os.path.exists(os.path.join(super_path, 'txt')):
    print('Creating txt folder...')
    os.mkdir(os.path.join(super_path, 'txt'))

if not os.path.isfile('../json/config.json'):
    print('Creating config.json...')
    while True:
        token = input('Enter your bot token: ')
        if token != '':
            break
        print('Please enter a valid token.')
    print('Saving config...')
    with open('../json/config.json', 'w') as f:
        json.dump({
            'token': token,
            'prefix': "!",
            'os': system
        }, f, indent=4)
    requirements = ["discord.py", "pubchempy", "rdkit", "pillow", "requests"]
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
    print('Done configuring...')


if __name__ == '__main__':
    if system == 'Linux':
        subprocess.run([venv_python, script_path], check=True)
    else:
        subprocess.run([sys.executable, script_path], check=True)
