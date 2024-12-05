""" new.py
The purpose of this program is to generate chemical structures locally instead of through Discord.
Copyright (C) 2024  github.com/brandongrahamcobb
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program.  
If not, see <https://www.gnu.org/licenses/>.
"""
import shlex
import traceback
import io
from PIL import Image

# A placeholder for your helpers module. Replace with your actual import statement.
import CobbBrandonGraham_helpers_051224 as helpers

def view_image(image: Image.Image, title: str):
    """Display the image using the default image viewer."""
    image.show(title=title)

def save_image(image: Image.Image, filename: str):
    """Save the image to a file."""
    image.save(filename)
    print(f"Image saved as {filename}")

def handle_option(option, molecules):
    try:
        if option == 'compare':
            if not molecules:
                print('No molecules provided.')
                return
            
            args = shlex.split(molecules)
            pairs = helpers.unique_pairs(args)  # Assuming this method generates unique pairs

            if not pairs:
                print('No valid pairs found.')
                return

            for pair in pairs:
                mol = helpers.get_mol(pair[0])
                refmol = helpers.get_mol(pair[1])

                if mol is None or refmol is None:
                    print(f'One or both of the molecules {pair[0]} or {pair[1]} are invalid.')
                    continue

                fingerprints = [
                    helpers.draw_fingerprint([mol, refmol]),
                    helpers.draw_fingerprint([refmol, mol])
                ]
                
                combined_image = helpers.combine(fingerprints, reversed(pair))

                # Check if combined_image is a BytesIO and convert it to an Image
                if isinstance(combined_image, io.BytesIO):
                    combined_image.seek(0)  # Move to the beginning of the BytesIO object
                    combined_image = Image.open(combined_image)  # Convert BytesIO to PIL Image
                
                # Save or show the image
                save_image(combined_image, 'molecule_comparison.png')
                view_image(combined_image, 'Molecule Comparison')

        elif option == 'glow':
            if not molecules:
                print('No molecules provided.')
                return
            
            args = shlex.split(molecules)
            fingerprints = []
            names = []
            molecule = helpers.get_mol(args[0])

            for _ in range(1):  # Assuming quantity is defined elsewhere
                names.append(args[0])
                fingerprints.append(helpers.draw_fingerprint([molecule, molecule]))

            combined_image = helpers.combine(fingerprints, names)

            # Check if combined_image is a BytesIO and convert it to an Image
            if isinstance(combined_image, io.BytesIO):
                combined_image.seek(0)  # Move to the beginning of the BytesIO object
                combined_image = Image.open(combined_image)  # Convert BytesIO to PIL Image

            # Save or show the image
            save_image(combined_image, 'molecule_glow.png')
            view_image(combined_image, 'Molecule Glow')

        elif option == 'shadow':
            if not molecules:
                print('No molecules provided.')
                return
            
            args = shlex.split(molecules)
            mol = helpers.get_mol(args[0])

            if mol is None:
                print('Invalid molecule name or structure.')
                return

            image = helpers.draw_watermarked_molecule(mol)  # Assuming this returns an Image object

            if isinstance(image, io.BytesIO):
                image.seek(0)  # Move to the beginning of the BytesIO object
                image = Image.open(image)  # Convert BytesIO to PIL Image

            # Save or show the shadow image
            save_image(image, f'{args[0]}.png')
            view_image(image, f'Shadow Molecule for {args[0]}')

        else:
            print('Invalid option. Use `compare`, `glow`, or `shadow`.')

    except Exception as e:
        print(f'An error occurred: {traceback.format_exc()}')

# Example usage
if __name__ == "__main__":
    while True:
        option = input("Enter option (compare, glow, shadow): ").strip()
        molecules = input("Enter molecules (space-separated): ").strip()
        handle_option(option, molecules)
