import pathlib
from functools import lru_cache
from multiprocessing import Lock, Process, Queue
from queue import Empty
from typing import Optional, Union

import pandas
from rdkit import Chem
from rdkit.Chem import GetSSSR, MolFromSmiles, MolToSmiles
from rdkit.Chem.rdchem import Mol as RDKitMol

# Create dataframe with hybridization list
ring_correction_path = (
    pathlib.Path(__file__).parent / "ring_corrections.csv"
)

def generate_hybridization_SMARTS(path_to_correction_CSV):
    """
    Adds ring SMARTS and corresponding RDKit molecules from the SMARTS to a ring correction CSV file, storing the information in a dataframe for future use.

    Parameters:
        path_to_correction_CSV (str): Path to the CSV file containing ring correction SMILES and their corresponding values.

    Returns:
        df (pandas.DataFrame): DataFrame containing information from the CSV, including the ring pattern SMARTS and the RDKit molecule generated using each SMARTS.
    """

    # Create molecule list
    mol_list = []

    # Read Ring Corrections into dataframe
    df = pandas.read_csv(path_to_correction_CSV)

    # Initialize list of values
    pattern_list = []
    pattern_mols = []

    # Loop through SMILES strings, create hybridization-specific SMARTS strings
    for smi in df['SMILES'].tolist():
        # Create mol, do NOT add explicit Hs, or there will be no match
        mol = Chem.MolFromSmiles(smi)
        mol_list.append(mol)
        
        # Create SMARTS pattern directly from the ring in the database
        pattern_string = Chem.MolToSmarts(mol)

        # Add Hs to molecule to get proper connectivity
        mol_with_Hs = Chem.AddHs(mol)

        # Create a connectivity list in order of atoms
        connectivity_list = []
        for _, atom in enumerate(mol_with_Hs.GetAtoms()):
            connectivity_list.append(len(atom.GetBonds()))

        # Use connectivity list to modify the pattern string based on number of brackets
        bracket_counter = 0
        updated_pattern_string = ""

        for char in pattern_string:
            if char == ']':
                updated_pattern_string += f';X{connectivity_list[bracket_counter]}]'
                bracket_counter += 1
            else:
                updated_pattern_string += char

        # # TESTING (for Shivani): Check if the pattern can be found in the source mol
        # match = mol.GetSubstructMatch(Chem.MolFromSmarts(updated_pattern_string))
        # try:
        #     assert(len(match) > 0)
        # except:
        #     print("ERROR")

        # Add to list
        pattern_list.append(updated_pattern_string)
        pattern_mols.append(Chem.MolFromSmarts(updated_pattern_string))

    # Set to column
    df['Hybridization'] = pattern_list
    df['SMARTS Mols'] = pattern_mols

    return df

# GLOBAL ring correction dataframe generated upon import
ring_correction_dataframe = generate_hybridization_SMARTS(ring_correction_path)

def _condense_rings(ringinfo_rings):
    """
    Identifies overlaps between sets in a list, particularly useful for creating ring index representations with multiple fused rings.

    Parameters:
        ringinfo_rings (tuple): Tuple of tuples containing unique rings, typically generated using RDKit [mol.GetRingInfo().AtomRings()].

    Returns:
        condensed_rings (list): List of individual rings represented as lists, with overlaps condensed.
    """

    # Type matching
    rings = list(ringinfo_rings)
    rings = [set(element) for element in rings]

    # Empty list to hold merged rings
    condensed_rings = []

    # Loop through all rings and condense
    for current_set in rings:
        new_set = current_set.copy()  # Create a copy to avoid modifying the original sets
        overlapping_indices = []

        for i, existing_set in enumerate(condensed_rings):
            if current_set.intersection(existing_set):
                overlapping_indices.append(i)

        for index in reversed(overlapping_indices):
            new_set |= condensed_rings.pop(index)

        condensed_rings.append(new_set)
    
    # Return condensed rings
    return [list(ele) for ele in condensed_rings]

def _get_direct_neighbors(mol, indices):
    """
    Returns a list of indices that includes the original indices along with their direct neighbors.

    Parameters:
        mol (RDKit mol): Larger molecule to which the indices belong.
        indices (list): Indices for which we want to find the direct neighbors.

    Returns:
        neighbors (list): List of original indices along with their direct neighbors.
    """

    # Set of neighbors (set to prevent overlap)
    neighbors = set()

    # Loop through neighbors of each index and add to list
    for idx in indices:
        atom = mol.GetAtomWithIdx(idx)

        for neighbor in atom.GetNeighbors():
            neighbors.add(neighbor.GetIdx())

    return list(neighbors)

def get_ring_correction_pathermo(mol):
    """
    Calculates the total ring correction for a given molecule.

    Parameters:
        mol (RDKit mol): Molecule of interest.
        ring_corr_dataframe (pandas.DataFrame): DataFrame containing ring correction values.

    Returns:
        total_ring_correction (float): Total combined ring correction value for all rings present in the molecule.
    """
    # Set flag
    missing_RC_flag = False

    # Extract relevant columns from df
    match_mols=ring_correction_dataframe['SMARTS Mols'].tolist()
    correction_values=ring_correction_dataframe['dHf_298'].tolist()

    # Get Ring Info for mol
    rings = _condense_rings(mol.GetRingInfo().AtomRings())

    # Intialize ring correction value
    total_ring_correction = 0

    # Loop through rings to find match
    for ring in rings:
        # Correction array for just the considered ring
        ring_corr = []

        # Add direct neighbors to list of main atoms in ring
        idx_list = _get_direct_neighbors(mol, ring)
        
        # Generate molecule from indices to match to database
        substructure_smiles = Chem.MolFragmentToSmiles(mol, idx_list, kekuleSmiles=True)
        substructure_mol = Chem.MolFromSmiles(substructure_smiles)

        # Search through all of the mols
        for i, m in enumerate(match_mols):
            match = substructure_mol.GetSubstructMatch(m)
            
            # If there is a match, add ring correction val
            if len(match) > 0:
                ring_corr.append(correction_values[i])
        
        # Currently adding the maximum, which is not always correct (need input from others here)
        if ring_corr: # Make sure ring exists in DB, might be taking max of nonexistent vals
            total_ring_correction += max(ring_corr)
        else:
            missing_RC_flag = True

    
    return round(total_ring_correction, 1), missing_RC_flag
    
