import json
import os
from typing import List, Dict

home = os.path.expanduser('~')
file = os.path.join(home, 'Downloads', "training.jsonl")

class TrainingFileBuilder:
    def __init__(self):
        """
        Initializes the TrainingFileBuilder with the specified training file.

        Args:
            training_file (str, optional): Path to the training JSONL file. Defaults to "training.jsonl".
        """
        self.training_file = file
        self.units = self._load_existing_units()

    def _load_existing_units(self) -> List[Dict[str, str]]:
        """
        Loads existing training units from the JSONL file.

        Returns:
            List[Dict[str, str]]: A list of existing units with 'prompt' and 'response'.
        """
        units = []
        if not os.path.exists(self.training_file):
            print(f"Training file '{self.training_file}' does not exist. A new file will be created.")
            return units

        try:
            with open(self.training_file, 'r', encoding='utf-8') as file:
                for line_num, line in enumerate(file, start=1):
                    line = line.strip()
                    if not line:
                        print(f"Skipping empty line #{line_num}.")
                        continue
                    try:
                        entry = json.loads(line)
                        prompt = entry.get('prompt', '').strip()
                        completion = entry.get('completion', '').strip()

                        if completion.startswith(' '):
                            completion = completion[1:]  # Remove leading space

                        if prompt and completion:
                            units.append({
                                "prompt": prompt,
                                "response": completion
                            })
                        else:
                            print(f"Skipping line #{line_num}: Missing 'prompt' or 'completion'.")
                    except json.JSONDecodeError as e:
                        print(f"Error decoding JSON on line #{line_num}: {e}")
        except IOError as e:
            print(f"An error occurred while reading the file: {e}")
        return units

    def add_responses(self, new_responses: List[str]) -> None:
        """
        Adds new assistant responses to the training units with unique prompt IDs.

        Args:
            new_responses (List[str]): A list of new assistant response strings.
        """
        next_id = self._get_next_prompt_id()
        for response in new_responses:
            if not isinstance(response, str):
                print(f"Skipping non-string response: {response}")
                continue
            response = response.strip()
            if not response:
                print("Skipping empty response.")
                continue
            unit = {
                "prompt": str(next_id),
                "response": response
            }
            self.units.append(unit)
            print(f"Added response with prompt ID {next_id}.")
            next_id += 1

    def _get_next_prompt_id(self) -> int:
        """
        Determines the next available prompt ID based on existing units.

        Returns:
            int: The next prompt ID.
        """
        if not self.units:
            return 1
        try:
            last_id = max(int(unit['prompt']) for unit in self.units)
            return last_id + 1
        except ValueError:
            print("Existing prompt IDs are not all integers. Starting prompt IDs from 1.")
            return 1

    def save(self) -> None:
        """
        Saves the current training units to the JSONL file.
        """
        try:
            with open(self.training_file, 'w', encoding='utf-8') as file:
                for unit in self.units:
                    json_line = json.dumps({
                        "prompt": unit["prompt"],
                        "completion": f" {unit['response']}"  # Add space before completion
                    }, ensure_ascii=False)
                    file.write(json_line + '\n')
            print(f"Training file '{self.training_file}' successfully updated with {len(self.units)} entries.")
        except IOError as e:
            print(f"An error occurred while writing to the file: {e}")
