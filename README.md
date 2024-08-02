# Lucy Bot

Lucy is a versatile Discord bot designed for managing molecular structures, interacting with emojis, and more. This README will guide you through setting up and running Lucy, including environment setup.

## Features

- **Colorize Role:** Change your role color with RGB values.
- **DMPurge:** Purge DM messages with the bot.
- **Draw:** Draw molecules or compare them graphically using RDKit.
- **Emoji Info:** Get information about emojis.
- **Load/Unload Cogs:** Dynamically load and unload bot cogs.
- **Purge:** Remove non-pinned messages from channels.
- **Reload:** Reload bot cogs.
- **SMILES:** Retrieve SMILES codes for chemicals.
- **Sync:** Synchronize the bot's command tree.

## Requirements

- Python 3.8+
- `discord.py` library
- `requests` library
- `asyncpraw` library
- `emoji` library
- `pubchempy` library
- `rdkit` library
- `pillow` library

## Installation

1. **Clone the repository:**

   ```bash
   git clone https://github.com/brandongrahamcobb/py.git
   cd py
   ```

2. **Run the setup script:**

   The `setup.py` script will handle environment setup, including creating a virtual environment, installing dependencies, and setting the bot token.

   ```bash
   python setup.py
   ```

   During setup, you will be prompted to enter your Discord bot token if it is not already present in `config.json`.

3. **Configure the bot:**

   - The setup script will create and update the `config.json` file in the `json` directory with your bot token.

4. **Create necessary directories:**

   The `setup.py` script will automatically create the required directories if they don't exist:

   - `../json/`
   - `../log/`
   - `../txt/`

## Running the Bot

To run the bot, execute the `main.py` script:

```bash
python main.py
```

This will start the bot and load the specified cogs.

## Code Explanation

- **Versioning:** The `get_version()` function manages version numbers and updates the version file.
- **Configuration Loading:** The `load_config()` function loads bot settings from `config.json`.
- **Logging:** The `setup_logging()` function configures logging to `../log/discord.log`.
- **Lucy Class:** A custom `commands.Bot` subclass that loads extensions and optionally syncs commands to a specific guild.
- **Main Function:** Initializes and starts the bot with the configured token and extensions.

## Commands

Here are the commands you can use with the bot:

- **!colorize <R> <G> <B>**: Change your role color with red, green, and blue values.
- **!dmpurge**: Purge DM messages with the bot (Luc).
- **!draw <MOLECULE> or !draw <MOLECULE1> <MOLECULE2>**: Draw a molecule or compare two molecules using RDKit.
- **!emoji <emoji>**: Get information about a given Unicode emoji character.
- **!load <extension>**: Load a cog.
- **!purge**: Delete non-pinned messages from the channel.
- **!reload <extension>**: Reload a cog.
- **!smiles <chemical>**: Retrieve the SMILES code for a chemical.
- **!unload <extension>**: Unload a cog.
- **!sync [~ | * | ^]**: Synchronize the bot's command tree globally or to the current guild.

## Environment Setup

The `setup.py` script performs the following:

- **System Update (Linux):** If running as root on Linux, updates the system using `pacman`.
- **Virtual Environment:** Creates a virtual environment if it doesn't exist.
- **Dependencies:** Installs required Python packages into the virtual environment.
- **Directory Creation:** Ensures necessary directories are created.
- **Token Setup:** Prompts for a Discord bot token if it is not present in `config.json`.

## Contributing

Contributions are welcome! Please fork the repository and submit a pull request with your changes.

## License

This project is licensed under the GPL-3.0 license. See the [LICENSE](LICENSE) file for details.

## Contact

For any questions or issues, please reach out to [brandongrahamcobb@icloud.com](mailto:brandongrahamcobb@icloud.com).
