from setuptools import setup, find_packages

setup(
    data_files=[
       ('py/bots/lucy', [
  #          'api/py-project.toml',
            'Roboto-Regular.ttf',
            'training.jsonl',
        ]),
        # Add other resource files or directories as needed
    ],
    entry_points={
        'console_scripts': [
            'lucy=bot.bot:main',
        ],
    },
    include_package_data=True,
    install_requires=[
        'asyncpg',
        'bs4',
        'discord.py',
        'emoji',
        'openai',
        'opencv-python',
        'pip',
        'pubchempy',
        'PyQt5',
        'requests',
        'rdkit',
        'selenium',
        'webdriver_manager',
    ],
    name='lucy',
    packages=find_packages(),
    python_requires='<3.13',
    version='4.8.6',
    zip_safe=False,
)
