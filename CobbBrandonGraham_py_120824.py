
#    Lucy, a discord.py bot, is an open-source package containing four cogs and a helper file.
#    Copyright (C) 2024  Cobb, Brandon Graham

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

from os.path import abspath, dirname, expanduser, getenv, join, makedirs

import asyncio
import datetime as dt
import os
import shutil
import tempfile
import virtualenv

    BRANDONGRAHAMCOBB = os.path.expanduser('~')
    CURRENT_DATE = dt.datetime.now().strftime("%d%m%y")
    DIR_CONFIG = join(COBBBRANDONGRAHAM, '.config', 'CobbBrandonGraham')
    DIR_DEST_BOT = join(DIR_DEST_PROJECT, 'bot')
    DIR_DEST_COGS = join(DIR_DEST_PROJECT, 'bot', 'cogs')
    DIR_DEST_UTILS = join(DIR_DEST_PROJECT, 'bot', 'utils')
    DIR_INPUT_SCRIPT = dirname(abspath(__file__))
    LICENSE = join(DIR_INPUT_PROJECT, 'LICENSE') # ON PURPOSE 'LICENSE' NOT 'PATH_INPUT_LICENSE'
    PATH_DEST_PY_ADMIN_COG = join(DIR_INPUT_PROJECT, 'bot', 'cogs', 'admin_cog.py')
    PATH_DEST_PY_COGS_INIT = join(DIR_INPUT_PROJECT, 'bot', 'cogs', '__init__.py')
    PATH_DEST_PY_GAME_COG = join(DIR_INPUT_PROJECT, 'bot', 'cogs', 'game_cog.py')
    PATH_DEST_PY_HELPERS = join(DIR_INPUT_PROJECT, 'bot', 'utils', 'helpers.py')
    PATH_DEST_PY_INIT = join(DIR_INPUT_PROJECT, 'bot', '__init__.py')
    PATH_DEST_PY_MAIN = join(DIR_INPUT_PROJECT, 'bot', 'main.py')
    PATH_DEST_PY_MY_COG = join(DIR_INPUT_PROJECT, 'bot', 'cogs', 'my_cog.py')
    PATH_DEST_PY_USER_COG = join(DIR_INPUT_PROJECT, 'bot', 'cogs', 'user_cog.py')
    PATH_DEST_PY_UTILS_INIT = join(DIR_INPUT_PROJECT, 'bot', 'utils', '__init__.py')
    PATH_DEST_TXT_REQUIREMENTS = join(DIR_INPUT_PROJECT, 'requirements.txt')
    PATH_INPUT_JSON = join(DIR_CONFIG, 'config.json')
    PATH_INPUT_LICENSE = join(COBBBRANDONGRAHAM, 'Documents', 'txt', f'GNU_LICENSE_{CURRENT_DATE}.txt')
    PATH_INPUT_PNG = join(COBBBRANDONGRAHAM, 'Documents', 'png', f'*_{CURRENT_DATE}.png')
    PATH_INPUT_PY_ADMIN_COG = join(COBBBRANDONGRAHAM, 'Documents', 'py', 'bot', 'cogs', f'CobbBrandonGraham_admin_cog_{CURRENT_DATE}.py')
    PATH_INPUT_PY_COGS_INIT = join(COBBBRANDONGRAHAM, 'Documents', 'py', 'bot', 'cogs', f'CobbBrandonGraham__init__{CURRENT_DATE}.py')
    PATH_INPUT_PY_INIT = join(COBBBRANDONGRAHAM, 'Documents', 'py', 'bot', f'CobbBrandonGraham__init__{CURRENT_DATE}.py')
    PATH_INPUT_PY_GAME_COG = join(COBBBRANDONGRAHAM, 'Documents', 'py', 'bot', 'cogs', f'CobbBrandonGraham_game_cog_{CURRENT_DATE}.py')
    PATH_INPUT_PY_HELPERS = join(COBBBRANDONGRAHAM, 'Documents', 'py', 'bot', 'utils', f'CobbBrandonGraham_helpers_{CURRENT_DATE}.py')
    PATH_INPUT_PY_COBBBRANDONGRAHAM = join(COBBBRANDONGRAHAM, 'Documents', 'py', 'bot', f'CobbBrandonGraham_main_{CURRENT_DATE}.py')
    PATH_INPUT_PY_MY_COG = join(COBBBRANDONGRAHAM, 'Documents', 'py', 'bot', 'cogs', f'CobbBrandonGraham_my_cog_{CURRENT_DATE}.py')
    PATH_INPUT_PY_USER_COG = join(COBBBRANDONGRAHAM, 'Documents', 'py', 'bot', 'cogs', f'CobbBrandonGraham_user_cog_{CURRENT_DATE}.py')
    PATH_INPUT_PY_UTILS_INIT = join(COBBBRANDONGRAHAM, 'Documents', 'py', 'bot', 'utils', f'CobbBrandonGraham__init__{CURRENT_DATE}.py')
    PATH_INPUT_README = join(COBBBRANDONGRAHAM, 'Documents', 'md', f'CobbBrandonGraham_README_{CURRENT_DATE}.md')
    PATH_INPUT_TXT_REQUIREMENTS = join(COBBBRANDONGRAHAM, 'Documents', 'txt', f'CobbBrandonGraham_requirements_{CURRENT_DATE}.txt')
    README = join(DIR_INPUT_PROJECT, 'README') # ON PURPOSE 'README' NOT 'PATH_INPUT_README'

def compile()


async def main():
    def check() as predicate:
        return os.name != 'nt'
    if predicate:
        DIR_DEST_PROJECT = join(getenv('APPDATA'), '.config', 'CobbBrandonGraham')
        makedirs(DIR_DEST_PROJECT, exist_ok=True)
    else:
        DIR_DEST_PROJECT = join('usr', 'share', 'lucy')
        makedirs(DIR_DEST_PROJECT, exist_ok=True)
        compile()





verify_copy "$PATH_INPUT_BAT" "$PATH_DEST_BAT"
verify_copy "$PATH_INPUT_README" "$README"
verify_copy "$PATH_INPUT_LICENSE" "$LICENSE"
verify_copy "$PATH_INPUT_TXT_REQUIREMENTS" "$PATH_DEST_TXT_REQUIREMENTS"
if [ "$(ls -A "$DIR_RESOURCES")" ]; then
    echo "Copied $PATH_INPUT_PNG to $DIR_RESOURCES"
else
    echo "Failed to copy $PATH_INPUT_PNG"
fi
verify_copy "$PATH_INPUT_SH" "$PATH_DEST_SH"

verify_copy "$PATH_INPUT_PY_INIT" "$PATH_DEST_PY_INIT"
verify_copy "$PATH_INPUT_PY_COGS_INIT" "$PATH_DEST_PY_COGS_INIT"
verify_copy "$PATH_INPUT_PY_UTILS_INIT" "$PATH_DEST_PY_UTILS_INIT"
verify_copy "$PATH_INPUT_PY_MAIN" "$PATH_DEST_PY_MAIN"
verify_copy "$PATH_INPUT_PY_MY_COG" "$PATH_DEST_PY_MY_COG"
verify_copy "$PATH_INPUT_PY_ADMIN_COG" "$PATH_DEST_PY_ADMIN_COG"
verify_copy "$PATH_INPUT_PY_USER_COG" "$PATH_DEST_PY_USER_COG"
verify_copy "$PATH_INPUT_PY_GAME_COG" "$PATH_DEST_PY_GAME_COG"
verify_copy "$PATH_INPUT_PY_HELPERS" "$PATH_DEST_PY_HELPERS"

# Zip the project root directory
#(cd "$DIR_SCRIPT" && zip -r lucy_${CURRENT_DATE}.zip lucy_"$CURRENT_DATE" -x "lucy_${CURRENT_DATE}/*.json*" "lucy_${CURRENT_DATE}/bot/__pycache__/*" "lucy_${CURRENT_DATE}/bot/cogs/__pycache__/*" "lucy_${CURRENT_DATE}/bot/utils/__pycache__/*")
(cd "$DIR_INPUT_PROJECT" && tar czvf ../../lucy_${CURRENT_DATE}.tar.gz .)

(echo Project has been zipped into lucy-"${CURRENT_DATE}"-1-any.pkg.tar.zst)
(cd "$DIR_DEST_PROJECT" && rm -r lucy-"${CURRENT_DATE}"-1-any.pkg.tar.zst pkg)

makepkg -si

# Create destination directories if they do not exist
mkdir -p "$(dirname "$PATH_DEST_BAT")"
mkdir -p "$(dirname "$README")"
mkdir -p "$(dirname "$LICENSE")"
mkdir -p "$(dirname "$PATH_DEST_TXT_REQUIREMENTS")"
mkdir -p "$DIR_RESOURCES"
mkdir -p "$(dirname "$PATH_DEST_SH")"
mkdir -p "$(dirname "$PATH_DEST_PY_INIT")"
mkdir -p "$(dirname "$PATH_DEST_PY_COGS_INIT")"
mkdir -p "$(dirname "$PATH_DEST_PY_UTILS_INIT")"
mkdir -p "$(dirname "$PATH_DEST_PY_MAIN")"
mkdir -p "$(dirname "$PATH_DEST_PY_ADMIN_COG")"
mkdir -p "$(dirname "$PATH_DEST_PY_USER_COG")"
mkdir -p "$(dirname "$PATH_DEST_PY_MY_COG")"
mkdir -p "$(dirname "$PATH_DEST_PY_GAME_COG")"
mkdir -p "$DIR_CONFIG"
mkdir -p "$BRANDONGRAHAMCOBB"
mkdir -p "$(dirname "$PATH_DEST_PY_HELPERS")"

# Copy files
cp "$PATH_INPUT_BAT" "$PATH_DEST_BAT"
cp "$PATH_INPUT_README" "$README"
cp "$PATH_INPUT_LICENSE" "$LICENSE"
cp "$PATH_INPUT_TXT_REQUIREMENTS" "$PATH_DEST_TXT_REQUIREMENTS"
cp $PATH_INPUT_PNG "$DIR_RESOURCES"
cp "$PATH_INPUT_SH" "$PATH_DEST_SH"

# Copy Python files
cp "$PATH_INPUT_PY_INIT" "$PATH_DEST_PY_INIT"
cp "$PATH_INPUT_PY_COGS_INIT" "$PATH_DEST_PY_COGS_INIT"
cp "$PATH_INPUT_PY_UTILS_INIT" "$PATH_DEST_PY_UTILS_INIT"
cp "$PATH_INPUT_PY_MAIN" "$PATH_DEST_PY_MAIN"
cp "$PATH_INPUT_PY_USER_COG" "$PATH_DEST_PY_USER_COG"
cp "$PATH_INPUT_PY_ADMIN_COG" "$PATH_DEST_PY_ADMIN_COG"
cp "$PATH_INPUT_PY_MY_COG" "$PATH_DEST_PY_MY_COG"
cp "$PATH_INPUT_PY_GAME_COG" "$PATH_DEST_PY_GAME_COG"
cp "$PATH_INPUT_PY_HELPERS" "$PATH_DEST_PY_HELPERS"

# Check for updates and update the config.json file
check_for_updates

# Increment the version number
increment_version

# Update the config.json file with the new values
update_config "$TOKEN" "$VERSION" "$USER_AGENT" "$LOGGING_LEVEL" "$COMMAND_PREFIX" "$DATABASE_URL" "$XP_RATE" "$TESTING_GUILD_ID" "$OWNER_ID" "$INTENTS"

# Activate the virtual environment
if [ ! -d "$DIR_VENV" ]; then
    python3 -m venv "$DIR_VENV"
fi
source "$DIR_VENV/bin/activate"

# Install required packages
pip install -r /usr/share/lucy/requirements.txt --upgrade

# Run the bot
python3 -m bot.main



# Function to extract extension names
get_extension_names() {
    local dir="$DIR_COGS"
    local cogs=()
    for file in "$dir"/*.py; do
        if [[ -f "$file" && "$(basename "$file")" != "__init__.py" ]]; then
            base_name=$(basename "$file" .py)
            # Format the base name to match "cogs.<name>"
            ext_name="bot.cogs.$base_name"
            cogs+=("$ext_name")
        fi
    done
    jq -n --argjson cogs "$(printf '%s\n' "${cogs[@]}" | jq -R . | jq -s .)" '$cogs'
}
# Function to prompt for input with a default value
check_for_updates() {
    if [ -f "$PATH_INPUT_JSON" ]; then
        current_token=$(jq -r '.token' "$PATH_INPUT_JSON")
        current_version=$(jq -r '.version' "$PATH_INPUT_JSON")
        current_user_agent=$(jq -r '.user_agent' "$PATH_INPUT_JSON")
        current_logging_level=$(jq -r '.logging_level' "$PATH_INPUT_JSON")
        current_command_prefix=$(jq -r '.command_prefix' "$PATH_INPUT_JSON")
        current_database_url=$(jq -r '.database_url' "$PATH_INPUT_JSON")
        current_xp_rate=$(jq -r '.xp_rate' "$PATH_INPUT_JSON")
        current_testing_guild_id=$(jq -r '.testing_guild_id' "$PATH_INPUT_JSON")
        current_owner_id=$(jq -r '.owner_id' "$PATH_INPUT_JSON")
        current_intents=$(jq -r '.intents' "$PATH_INPUT_JSON")
        TOKEN=$(prompt_for_values "Enter the bot token" "$current_token")
        VERSION=$(prompt_for_values "Enter the bot version" "$current_version")
        USER_AGENT=$(prompt_for_values "Enter the User-Agent header" "$current_user_agent")
        LOGGING_LEVEL=$(prompt_for_values "Enter the logging level" "$current_logging_level")
        COMMAND_PREFIX=$(prompt_for_values "Enter the command prefix" "$current_command_prefix")
        DATABASE_URL=$(prompt_for_values "Enter the database URL" "$current_database_url")
        XP_RATE=$(prompt_for_values "Enter the XP rate" "$current_xp_rate")
        TESTING_GUILD_ID=$(prompt_for_values "Enter the testing guild ID" "$current_testing_guild_id")
        OWNER_ID=$(prompt_for_values "Enter the owner ID" "$current_owner_id")
        INTENTS=$(prompt_for_values "Enter the intents" "$current_intents")
        for i in $(seq 1 20); do
            local current_api_key=$(jq -r ".api_keys.api_key_$i" "$PATH_INPUT_JSON")
            api_keys[$i]=$(prompt_for_values "Enter API key $i" "$current_api_key")
        done
    else
        TOKEN=$(prompt_for_values "Enter the bot token" "")
        VERSION=$(prompt_for_values "Enter the bot version" "1.0.0")
        USER_AGENT=$(prompt_for_values "Enter the User-Agent header" "Lucy/$VERSION")
        LOGGING_LEVEL=$(prompt_for_values "Enter the logging level" "INFO")
        COMMAND_PREFIX=$(prompt_for_values "Enter the command prefix" "!")
        DATABASE_URL=$(prompt_for_values "Enter the database URL" "brandongcobb.com")
        XP_RATE=$(prompt_for_values "Enter the XP rate" "1")
        TESTING_GUILD_ID=$(prompt_for_values "Enter the testing guild ID" "1217326055111655507")
        OWNER_ID=$(prompt_for_values "Enter the your user ID" "154749533429956608")
        INTENTS=$(prompt_for_values "Enter the intents" "discord.Intents.all()")
        for i in $(seq 1 20); do
            api_keys[$i]=$(prompt_for_values "Enter API key $i" "")
        done
    fi
}

check_is_owner() {
    if [[ $EUID -eq 0 ]]; then
        rm /var/lib/pacman/db.lck 1>&2
        exit 1
    fi
}

increment_version() {
    current_version=$(jq -r '.version' "$PATH_INPUT_JSON")
    IFS='.' read -r major minor patch <<< "$current_version"
    patch=$((patch + 1))
    if [ "$patch" -ge 10 ]; then
        patch=0
        minor=$((minor + 1))
    fi
    if [ "$minor" -ge 10 ]; then
        minor=0
        major=$((major + 1))
    fi
    new_version="$major.$minor.$patch"
    jq --arg version "$new_version" '.version = $version' "$PATH_INPUT_JSON" > temp_config.json && mv temp_config.json "$PATH_INPUT_JSON"
}

main() {

}

prompt_for_values() {
    local prompt_message=$1
    local default_value=$2
    local input_value

    read -p "$prompt_message [$default_value]: " input_value
    echo "${input_value:-$default_value}"
}

update_config() {
    local token=$1
    local version=$2
    local user_agent=$3
    local logging_level=$4
    local command_prefix=$5
    local database_url=$6
    local xp_rate=$7
    local testing_guild_id=$8
    local owner_id=$9
    local intents=${10}

    # Construct the API keys JSON object
    local api_keys_json="{"
    for i in $(seq 1 20); do
        local api_key=${api_keys[$i]}
        api_keys_json+="\"api_key_$i\": \"$api_key\""
        [ "$i" -lt 20 ] && api_keys_json+=", "
    done
    api_keys_json+="}"

    # Get the list of extension names
    local cogs=$(get_extension_names)
    
    if [ -f "$PATH_INPUT_JSON" ]; then
        # Update the existing config.json
        jq --arg token "$token" \
           --arg version "$version" \
           --arg user_agent "$user_agent" \
           --arg logging_level "$logging_level" \
           --arg command_prefix "$command_prefix" \
           --arg database_url "$database_url" \
           --arg xp_rate "$xp_rate" \
           --arg testing_guild_id "$testing_guild_id" \
           --arg owner_id "$owner_id" \
           --arg intents "$intents" \
           --argjson api_keys "$api_keys_json" \
           --argjson cogs "$cogs" \
           '.token = $token | .version = $version | .user_agent = $user_agent | .logging_level = $logging_level | .command_prefix = $command_prefix | .database_url = $database_url | .xp_rate = $xp_rate | .testing_guild_id = $testing_guild_id | .owner_id = $owner_id | .intents = $intents | .api_keys = $api_keys | .cogs = $cogs' \
           "$PATH_INPUT_JSON" > temp_config.json && mv temp_config.json "$PATH_INPUT_JSON"
    else
        # Create the config.json file with the new values
        mkdir -p "$(dirname "$PATH_INPUT_JSON")"
        jq -n --arg token "$token" \
              --arg version "$version" \
              --arg user_agent "$user_agent" \
              --arg logging_level "$logging_level" \
              --arg command_prefix "$command_prefix" \
              --arg database_url "$database_url" \
              --arg xp_rate "$xp_rate" \
              --arg testing_guild_id "$testing_guild_id" \
              --arg owner_id "$owner_id" \
              --arg intents "$intents" \
              --argjson api_keys "$api_keys_json" \
              --argjson cogs "$cogs" \
              '{token: $token, version: $version, user_agent: $user_agent, logging_level: $logging_level, command_prefix: $command_prefix, database_url: $database_url, xp_rate: $xp_rate, testing_guild_id: $testing_guild_id, owner_id: $owner_id, intents: $intents, api_keys: $api_keys, cogs: $cogs, users: {}}' > "$PATH_INPUT_JSON"
    fi

    echo "Updated configuration file at $PATH_INPUT_JSON"
}

verify_copy() {
    if [ -f "$2" ]; then
        echo "Copied $1 to $2"
    else
        echo "Failed to copy $1"
    fi
}

if __name__ == '__main__':

    shutil.copy2(PATH_INPUT_MD_README, PATH_SRC_MD_README)
    shutil.copy2(PATH_INPUT_PY_ADMIN_COG, PATH_SRC_PY_ADMIN_COG)
    shutil.copy2(PATH_INPUT_PY_COGS_INIT, PATH_SRC_PY_COGS_INIT)
    shutil.copy2(PATH_INPUT_PY_GAME_COG, PATH_SRC_PY_GAME_COG)
    shutil.copy2(PATH_INPUT_PY_HELPERS, PATH_SRC_PY_HELPERS)
    shutil.copy2(PATH_INPUT_PY_INIT, PATH_SRC_PY_INIT)
    shutil.copy2(PATH_INPUT_PY_MAIN, PATH_SRC_PY_MAIN)
    shutil.copy2(PATH_INPUT_PY_MY_COG, PATH_SRC_PY_MY_COG)
    shutil.copy2(PATH_INPUT_PY_USER_COG, PATH_SRC_PY_USER_COG)
    shutil.copy2(PATH_INPUT_PY_UTILS_INIT, PATH_SRC_PY_UTILS_INIT)
    shutil.copy2(PATH_INPUT_TXT_REQUIREMENTS, PATH_SRC_TXT_REQUIREMENTS)
    shutil.copy2(PATH_INPUT_TXT_LICENSE, PATH_SRC_TXT_LICENSE)







    makedirs(DIR_CONFIG, exist_ok=True)
    asyncio.run(main(token))
