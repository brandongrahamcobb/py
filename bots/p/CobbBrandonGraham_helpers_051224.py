''' helpers.py
    Copyright (C) 2024 github.com/brandongrahamcobb

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''


from datetime import datetime as dt
from io import BytesIO
from os import makedirs
from os.path import abspath, dirname, exists, expanduser, isfile, join
from PIL import Image, ImageFont, ImageDraw
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, DataStructs, Draw, rdDepictor, rdFingerprintGenerator, rdFMCS
#rdDepictor.SetPreferCoordGen(True)
from rdkit.Chem.Draw import rdMolDraw2D, SimilarityMaps
from rdkit.DataStructs import FingerprintSimilarity, TanimotoSimilarity
from typing import List

import colorsys
import itertools
import math
import os
import pubchempy as pcp
import traceback

current_date = dt.now().strftime('%d%m%y')
dir_base = abspath(__file__)

path_home = expanduser('~')
path_base = dirname(dir_base)

def add_watermark(image: BytesIO, watermark_text: str = '~spooky~') -> BytesIO:
    RGB_image = Image.open(image)
    RGBA_image = RGB_image.convert('RGBA')
    draw = ImageDraw.Draw(RGBA_image)
    width, height = RGBA_image.size
    diagonal = math.sqrt(width**2 + height**2)
    font_size = int(diagonal / 15)
    try:
        font = ImageFont.truetype('Roboto-Regular.ttf', font_size)  # Replace with the path to your font file
    except IOError:
        font = ImageFont.load_default()
    while True:
        bbox = draw.textbbox((0, 0), watermark_text, font=font)
        text_width = bbox[2] - bbox[0]
        if text_width <= 512 or font_size <= 1:
            break
        font_size -= 1
        font = ImageFont.truetype('Roboto-Regular.ttf', font_size)  # Replace with the path to your font file
    text_height = bbox[3] - bbox[1]
    text_x = (width - text_width) / 2
    text_y = height - (2 * text_height)
    watermark_image = Image.new('RGBA', RGBA_image.size, (0, 0, 0, 0))
    watermark_draw = ImageDraw.Draw(watermark_image)
    watermark_draw.text((text_x, text_y), watermark_text, font=font, fill=(255, 255, 255, 64))
    mask = watermark_image.split()[3]
    RGBA_image.paste(watermark_image, (0, 0), mask)
    output = BytesIO()
    RGBA_image.save(output, format= 'PNG')
    output.seek(0)
    return output

def adjust_hue_and_saturation(image, hue_shift, saturation_shift) -> BytesIO:
    image = image.convert('RGB')
    pixels = list(image.getdata())
    adjusted_pixels = []
    for r, g, b in pixels:
        h, s, v = colorsys.rgb_to_hsv(r / 255.0, g / 255.0, b / 255.0)
        h = (h + hue_shift / 360.0) % 1.0
        s = min(max(s + saturation_shift / 100.0, 0), 1)
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        adjusted_pixels.append((int(r * 255), int(g * 255), int(b * 255)))
    new_image = Image.new('RGB', image.size)
    new_image.putdata(adjusted_pixels)
    output = BytesIO()
    new_image.save(output, format= 'PNG')
    output.seek(0)
    return output

def combine(images: list, names: list) -> BytesIO:
    combined_images = []
    for index, (bytes_io, name) in enumerate(zip(images, names)):
        img = Image.open(bytes_io)
        inverted_image = Image.eval(img, lambda x: 255 - x)
        img_bytes = BytesIO()
        inverted_image.save(img_bytes, format='PNG')
        img_bytes_final = add_watermark(img_bytes, watermark_text=name)
        img_final = Image.open(img_bytes_final)
        combined_images.append(img_final)
    widths, heights = zip(*(img.size for img in combined_images))
    total_width = sum(widths)
    max_height = max(heights)
    combined_img = Image.new('RGB', (total_width, max_height))
    x_offset = 0
    for img in combined_images:
        combined_img.paste(img, (x_offset, 0))
        x_offset += img.width
    output = adjust_hue_and_saturation(combined_img, hue_shift=-180, saturation_shift=160)
    output.seek(0)
    return output

def draw_fingerprint(pair) -> BytesIO:
#    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, countSimulation=True)
#    def get_fp(mol, *args, **kwargs):
#        return mfpgen.GetFingerprint(mol)
    d2d = rdMolDraw2D.MolDraw2DCairo(1024, 1024)
    d2d.prepareMolsBeforeDrawing = False
    Options = d2d.drawOptions()
    Options.prepareMolsBeforeDrawing = False
    Options.includeMetadata = False
    Options.bondLineWidth = 4.0
    d2d.SetDrawOptions(Options)
    mol1 = rdMolDraw2D.PrepareMolForDrawing(pair[0], kekulize=True)
    mol1.UpdatePropertyCache(False)
    mol2 = rdMolDraw2D.PrepareMolForDrawing(pair[1], kekulize=True)
    mol2.UpdatePropertyCache(False)
    fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(mol1, mol2, lambda m, i: SimilarityMaps.GetMorganFingerprint(m, i, radius=2, fpType='bv', nBits=8192), draw2d=d2d, drawingOptions=Options)
#    fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(pair[1], pair[0], get_fp, draw2d=d2d) #colorMap=brighter_color, draw2d=d2d)
    d2d.FinishDrawing()
    drawing = d2d.GetDrawingText()
    output = BytesIO(drawing)
    return output

def draw_watermarked_molecule(molecule) -> BytesIO:
    resolved_name = get_molecule_name(molecule)
    d2d = rdMolDraw2D.MolDraw2DCairo(128, 128)
    Options = d2d.drawOptions()
    Options.prepareMolsBeforeDrawing = False
    Options.includeMetadata = False
    Options.bondLineWidth = 4.0
    d2d.SetDrawOptions(Options)
    rdMolDraw2D.SetDarkMode(Options)
    mol = rdMolDraw2D.PrepareMolForDrawing(molecule, kekulize=True)
    mol.UpdatePropertyCache(False)
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    drawing = d2d.GetDrawingText()
    output = add_watermark(BytesIO(drawing), watermark_text=resolved_name)
    return output

def get_proximity(default, input) -> float:
    default_fp = AllChem.GetMorganFingerprintAsBitVect(default, 2)
    input_fp = AllChem.GetMorganFingerprintAsBitVect(input, 2)
    similarity = DataStructs.FingerprintSimilarity(default_fp, input_fp)
    return similarity

def get_mol(arg):
    try:
        compounds = pcp.get_compounds(arg, 'name')
        compound_data = compounds[0].to_dict(properties=['isomeric_smiles'])
        mol = Chem.MolFromSmiles(compound_data['isomeric_smiles'])
        if mol is None:
            raise ValueError('Invalid SMILES string')
        return mol
    except:
        mol = Chem.MolFromSmiles(arg)
        if mol is None:
            raise ValueError('Invalid SMILES string')
        return mol

def get_molecule_name(molecule) -> str:
    smiles = Chem.MolToSmiles(molecule)
    if smiles:
        compounds = pcp.get_compounds(smiles, 'smiles') #, record_type='3d')
        compound_data = compounds[0].to_dict(properties=['synonyms'])
        if not compounds:
            raise ValueError('No compound found for the given SMILES string')
        return compound_data['synonyms'][0]
    return 'Unknown'

def unique_pairs(strings_list):
    pairs = list(itertools.combinations(strings_list, 2))
    sorted_pairs = [sorted(list(pair)) for pair in pairs]
    sorted_pairs_overall = sorted(sorted_pairs)
    return sorted_pairs_overall
