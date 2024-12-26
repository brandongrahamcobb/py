''' get_mol.py  The purpose of this program is fetch a molecule SMILES from pubchem from a name from cd ../.
    Copyright (C) 2024  github.com/brandongrahamcobb

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
from rdkit import Chem

import pubchempy as pcp

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
