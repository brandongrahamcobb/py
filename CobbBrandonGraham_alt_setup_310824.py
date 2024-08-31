from setuptools import find_packages, setup

setup(
    name='custom_bot',
    version='1.0.0',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'discord.py',
        'gradio_client',
        'pip',
    ],
    entry_points={
        'console_scripts': [
            'custom_bot=bot.main:run',
        ],
    },
    zip_safe=False,
)
