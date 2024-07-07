import subprocess
import venv
import json
import sys
import os

location = os.path.dirname(os.path.abspath(__file__))
system = 'Windows' if os.name == 'nt' else 'Linux'
venv_path = os.path.join(location, 'venv') if system == 'Linux' else None
venv_python = os.path.join(venv_path, 'bin', 'python') if system == 'Linux' else None
script_path = os.path.join(location, 'Main.py')

print(f'System: {system}')
print(f'Location: {location}')

def clear():
    if system == 'Windows': os.system('cls')
    elif system == 'Linux': os.system('clear')        

if not os.path.isfile('config.json'):
    print('Creating config.json...')
    
    while True:
        token = input('Enter your bot token: ')
        if token != '':
            break
        print('Please enter a valid token.')
    
    print('Saving config...')
    
    with open('config.json', 'w') as f:
        json.dump({
            'token': token,
            'prefix': "!",
            'os': system
        }, f, indent=4)
        
    # add dependencies, im not sure which ones you've had before
    requirements = ["discord.py", "pubchempy", "rdkit", "pillow", "requests"]
        
    if system == 'Linux':
        print(f'Python Location: {sys.executable}')
        print(f'Python Version: {sys.version}')
        q = input('Do you want to create a venv? (y/n): ')
        if q.lower() == 'y':
            print('Creating virtual environment...')
            venv.create(venv_path, with_pip=True)
            print('Activating virtual environment...')
            subprocess.run([venv_python, '-m', 'venv', 'activate'], check=True)
            print('Installing dependencies...')
            subprocess.run([venv_python, '-m', 'pip', 'install'] + requirements, check=True)
            print('Virtual environment created successfully.')
    else:
        print('Installing dependencies...')
        subprocess.run([sys.executable, '-m', 'pip', 'install'] + requirements, check=True)
    
    print('Done configuring...')


if __name__ == '__main__':
    subprocess.run([venv_python, script_path], check=True)
