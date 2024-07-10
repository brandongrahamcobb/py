
# ChemBot

![ChemBot Logo](path/to/logo.png)

**ChemBot** is a powerful Discord bot designed to assist chemists with various tasks, including molecular visualization, chemical property prediction, and much more. By leveraging the power of RDKit and discord.py, ChemBot brings advanced chemistry tools directly to your Discord server.

## Features

- **Molecular Visualization**: Generate 2D and 3D visualizations of molecules.
- **Property Prediction**: Calculate molecular properties such as molecular weight, LogP, and more.
- **Chemical Structure Search**: Search for molecules by structure or substructure.
- **Reaction Prediction**: Predict chemical reactions based on input molecules.
- **Database Integration**: Integrate with various chemical databases for extended functionality.
- **User-Friendly Commands**: Easy-to-use commands with detailed help options.

## Installation

### Prerequisites

- Python 3.7 or higher
- RDKit
- discord.py

### Setup

1. **Clone the repository**
   ```bash
   git clone https://github.com/yourusername/ChemBot.git
   cd ChemBot
   ```

2. **Create a virtual environment**
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows, use `venv\Scripts\activate`
   ```

3. **Install dependencies**
   ```bash
   pip install -r requirements.txt
   ```

4. **Set up your Discord bot token**
   - Create a `.env` file in the project root directory:
     ```
     DISCORD_TOKEN=your_discord_bot_token
     ```

5. **Run the bot**
   ```bash
   python bot.py
   ```

## Usage

Invite ChemBot to your server using [this link](https://discord.com/oauth2/authorize?client_id=your_client_id&scope=bot&permissions=your_permissions).

### Commands

- `!mol <SMILES>`: Generate a 2D visualization of a molecule from its SMILES string.
- `!props <SMILES>`: Calculate and display properties of a molecule.
- `!help`: Display the list of available commands.

### Examples

- Generate a 2D visualization of benzene:
  ```
  !mol c1ccccc1
  ```

- Calculate properties of ethanol:
  ```
  !props CC(O)C
  ```

## Contributing

We welcome contributions! Please read our [Contributing Guide](CONTRIBUTING.md) to get started.

## License

ChemBot is licensed under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for more details.

## Acknowledgements

- [RDKit](https://www.rdkit.org/)
- [discord.py](https://discordpy.readthedocs.io/en/stable/)

![ChemBot Demo](path/to/demo.gif)
