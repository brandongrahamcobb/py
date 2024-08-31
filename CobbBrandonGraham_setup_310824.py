from setuptools import find_packages, setup

setup(
    name='lucy',
    version='300824',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'biopython',
        'bs4',
        'discord.py',
        'emoji',
        'google-api-python-client',
        'googletrans==4.0.0-rc1',
        'pip',
        'pubchempy',
        'pillow',
        'requests',
        'rdkit',
        'selenium',
        'webdriver_manager',
    ],
    entry_points={
        'console_scripts': [
            'lucy=bot.main:run',
        ],
    },
    data_files=[
        ('/usr/share/lucy/resources', [
            'resources/config.json',
            'resources/users.json',
        ]),
        # Add other resource files or directories as needed
    ],
    zip_safe=False,
)
