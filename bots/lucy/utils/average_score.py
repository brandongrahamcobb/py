import asyncio
import json
import os

async def main():

    home = os.path.expanduser('~') 
 
    path_training = os.path.join(home, 'py', 'bots', 'lucy', 'training.jsonl')

    with open(path_training, 'r') as file:
        taining_data = json.load(file)
        print(training_data)

    file.close()

if __name__ == '__main__':

    asyncio.run(main())
