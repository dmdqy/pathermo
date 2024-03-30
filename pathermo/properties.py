from rdkit.Chem import MolFromSmiles, MolFromSmarts, BondType, CanonSmiles
from functools import lru_cache
from pathlib import Path
import pandas as pd
import json
import time
from .ring_functions import get_ring_correction_pathermo

group_values_path = Path(__file__).parent / "group_values"
molecules_path = Path(__file__).parent / "molecules"
ring_correction_path = Path(__file__).parent / "ring_corrections.csv"

with open(group_values_path, encoding="utf8") as f:                          # file with group values
        lines = f.readlines()          
group_Hf = dict()                                                              # {"C 4 C H H H": float}        
for line in lines[1:]:                                                         # skip first line
    if line[0] != "/" and line[0:2] != '\n' and line[0] != ' ':                # skip comments and empty lines
        if line[0:3] == "END":                                                 # stop at END
            break        
        parts = line.split() 
        try:
            group_Hf[" ".join(parts[0 : int(parts[1]) + 2])] = float(parts[int(parts[1]) + 2])  # add key and value to dict
        except:
            group_Hf[" ".join(parts[0 : int(parts[1]) + 2])] = float(parts[int(parts[1]) + 2][0:-1]) # if value has one unexpected character at the end                            

with open(molecules_path) as g:                                          # molecule file, stores Hf of simple molecules
        lines = g.readlines()          
molecules_Hf = dict()                                                          # {"O": float} Hf data        
for line in lines:                                                             # skip first line
    if line[0] != "/" and line[0:2] != '\n' and line[0] != ' ':                # skip comments and empty lines      
        parts = line.split()         
        molecules_Hf[parts[0]] = float(parts[1])                               # add key and value to dict


@lru_cache(16384) # cache, 2^14 = 16384
def Hf(smiles, return_missing_groups = False, return_regardless_of_missing_RCs = True):    
           
    if smiles in molecules_Hf:                    
        return molecules_Hf[smiles]
    
    canon_smiles = CanonSmiles(smiles)    
    if canon_smiles in molecules_Hf:
        return molecules_Hf[canon_smiles]
             
    polyatomic_groups = {                                                      # atom neighboring terminal heavy atoms are treated as a whole
        "[CX3+0]=[O+0]": "CO",     # *-C(=O)-*      start with the base atom where it connects to other groups
        "[Cv3+0]=[O+0]": "CO.",     # *-C(=O)  radical
        "[SX3+0]=[O+0]": "SO",     # *-S(=O)-*
        "[SX4+0](=[O+0])=[O+0]": "So",      # *-S(=O)(=O)-*
        "[NX3+1](=[O+0])[O-1H0]": "No",     # *-NO2
        "[NX2+0]=[O+0]": "NO",      # *-N=O
        "[CX3+0]=[SX1+0]": "CS",   # *-C(=S)-*
        "[PX4+0]=[O+0]": "PO",     # *-P(=O)(-*)-*   
        "[CX2+0]#[N+0]": "CN",     # *-C#N
        # PN                       # lacks data
        #Bo=BO3                    # lacks data
        "[NX2+0]=[C+0]=[O+0]": "Nco",     # *-N=C=O, isocyanate group
        "[C+0](=[NX2+0])=[O+0]": "Nco",   # same group as above, write middle O as base so it's treated as part of the dependent group.
        }
    
    polyatomic_groups_but_dependent = {    # some groups are not independent, they appear together with a carbon atom like C 4 H H H No
        "No",      # -NO2
        "NO",      # -N=O
        "CN",      # -C#N
        "CO.",     # *-C(=O) radical
        "Nco",      # *-N=C=O, isocyanate group
        }
    
    m = MolFromSmiles(smiles)                                             # rdkit mol
    polyatomic_groups_idx = dict()                                             # {int: [str, int]}  atom idx : [group name, terminal atom idx]
    for SMARTS in polyatomic_groups:
        matches = m.GetSubstructMatches(MolFromSmarts(SMARTS))            # find if such groups are in mol
        for tup in matches:        
            polyatomic_groups_idx[tup[0]] = [polyatomic_groups[SMARTS]] + list(tup[1:]) # put key and value to dict
    
    total_Hf = 0                                                               # Hf of molecule    
    double_bond_set = {element for tupl in m.GetSubstructMatches(MolFromSmarts("[*]=[*]")) for element in tupl} # contains all atom with double bond, CD for C=C, N=C,..     
    triple_bond_set = {element for tupl in m.GetSubstructMatches(MolFromSmarts("[*]#[*]")) for element in tupl} # contains all atom with triple bond, CT for C#C       
    fused_aromatic = {element for tupl in m.GetSubstructMatches(MolFromSmarts("[a$(*(:a)(:a):a)]")) for element in tupl} # contains all fused aromatic atoms, CF in groups file, or CBF in Benson book

    group_checked = False                                                      # store if at least one group is used

    missing_groups = list()
    missing_groups_flag = False

    for atom in m.GetAtoms():                                                  # iterate over all atoms in mol  

        total_degree = atom.GetTotalDegree()                                   # number of neighbors, including H
        if total_degree == 1:                                                  # skip terminal atom 
            continue       
        atom_symbol = atom.GetSymbol()                                         # get atom symbol, C, O, N,...
        poly_flag = False                                                      # mark if this atom is the base of a polyatomic group
        atom_idx = atom.GetIdx()                                               # idx of this atom
        if atom_idx in polyatomic_groups_idx:
            atom_symbol = polyatomic_groups_idx[atom_idx][0]                   # if atom is base of a group, assign the group name to it
            if atom_symbol in polyatomic_groups_but_dependent:                 # if this group is not independent, skip
                continue
            total_degree = total_degree - len(polyatomic_groups_idx[atom_idx][1:]) # get degree of group
            poly_flag = True
        elif atom.GetFormalCharge() != 0:
            atom_symbol = atom_symbol + "+"                                    # assign charge
        elif atom.GetNumRadicalElectrons() != 0:                               # check if radical
            atom_symbol = atom_symbol + "." 
                
        neighbors = atom.GetTotalNumHs() * ["H"]                               # list of atom's neighors, initiate with H
        aromatic_flag = atom.GetIsAromatic()                                   # mark if atom is aromatic
        unsaturated = 0                                                        # to use if atom has double/triple bonds
        
        for x in atom.GetNeighbors():                                          # iterate over atom's neighbors
            x_idx = x.GetIdx()                                                 # idx of this neighbor
            # next add atom's neighbors to neighbors list
            if not poly_flag or x_idx not in polyatomic_groups_idx[atom_idx][1:]: # skip x if atom is in a poly group and x is a terminal
                if x_idx in polyatomic_groups_idx:
                    neighbors.append(polyatomic_groups_idx[x_idx][0])          # if x is the base of a poly group, add group name as neighbor                         
                else:     
                    x_symbol = x.GetSymbol()                                   # symbol of x               
                    if m.GetBondBetweenAtoms(atom_idx, x_idx).GetBondType() == BondType.SINGLE: # if bond between atom and x is single
                        
                        if x.GetFormalCharge() != 0:
                            x_symbol = x_symbol + "+"                          # assign charge                       
                        if x.GetNumRadicalElectrons() != 0:                    # check if radical
                            x_symbol = x_symbol + "." 
                        
                        if x.GetIsAromatic():
                            neighbors.append(x_symbol + "B")                   # if x is aromtic, add xB as its name
                        else:
                            if x_idx in double_bond_set:
                                if x_symbol != "N":
                                    if x_symbol[-1] != ".":
                                        neighbors.append(x_symbol + "D")           # if x has double bonds to others, add xD as its name
                                    else:
                                        neighbors.append(x_symbol[:-1] + "D.")
                                else:
                                    neighbors.append(x_symbol + "I")           # imino N
                            elif x_idx in triple_bond_set:
                                neighbors.append(x_symbol + "T")               # if x has triple bond to others, add xT as its name                           
                            else:
                                neighbors.append(x_symbol)                     # normal condition, add x's symbol as its name
                    elif m.GetBondBetweenAtoms(atom_idx, x_idx).GetBondType() == BondType.DOUBLE: # if double bond between atom and x
                        unsaturated += 1                                       # x don't count as a neighbor, degree -1
                        if atom_symbol != "N":
                            if atom_symbol[-1] != ".":
                                atom_symbol += "D"                                 # add D to atom name, example: CD for C in C=C
                            else:
                                atom_symbol = atom_symbol[:-1] + "D."
                            if atom_symbol == "CDD":
                                atom_symbol = "CA"                             # allenic C, middle C of C=C=C, rare case
                        else:                                                  # imino, azo N cases
                            if x_symbol == "C":                                # imino
                                atom_symbol = "NI"
                            elif x_symbol == "N":                              # azo
                                atom_symbol = "NA"
                                                        
                    elif m.GetBondBetweenAtoms(atom_idx, x_idx).GetBondType() == BondType.TRIPLE: # if triple bond between atom and x
                        unsaturated += 1                                       # x don't count as a neighbor, degree -1
                        atom_symbol += "T"                                     # add T to atom name, example: CT for C in C#C
                                            
        # next deal with situation when atom is aromatic        
        if aromatic_flag and atom_idx not in fused_aromatic:                   # if atom is aromatic but not fused atom
            if atom_symbol == "C":
                total_degree = 1                                               # aromatic C count only one neighbor, except for special cases
                atom_symbol = "CB"                                             # assign atom name as CB                       
            
            elif atom_symbol == "O" or atom_symbol == "S":                     # if aromatic O or S in a furan ring            
                neighbors = ["CB", "CB"]                                       # Could be N but ignore

            elif atom_symbol == "N":                                           # if aromatic N    
                atom_symbol = "NB"                                             # assign atom name as NB  
                total_degree = total_degree - 2                                # aromatic N may have one neighbor for pyrrole or 0 for pyridine

            elif atom_symbol == "CD":                                          # rare situations like O=C1C=COC=C1, could be c=O, c=C, c=N, c=S
                for xx in atom.GetNeighbors():
                    if not xx.GetIsAromatic():                                 # deal with the non-aromatic neighbor
                        xx_symbol = xx.GetSymbol()
                        if xx_symbol == "O":
                            atom_symbol = "CO"
                            neighbors = ["CB", "CB"]                           # could be others but ignore

            elif atom_symbol == "C.":                                          # rare cases for aromatic radicals
                total_degree = 0                                               #
                atom_symbol = "CB."                                            #
                
         
        elif atom_idx in fused_aromatic:                                       # if atom is fused aromatic 
            atom_symbol += "F"                                                 # name CF
            for xx in atom.GetNeighbors():
                if xx.GetIdx() in fused_aromatic:                              # check if neighbor is also fused C
                    neighbors.append("CF")
                else:
                    if xx.GetSymbol() == "C":
                        neighbors.append("CB")                                 # could be O but ignore
                    elif xx.GetSymbol() == "N":
                        neighbors.append("NB")
                    elif xx.GetSymbol() == "O":
                        neighbors.append("O")
  
        neighbors.sort()                                                       # sort neighbors to match with group data
        key = " ".join([atom_symbol, str(total_degree - unsaturated)] + neighbors) # construct key for looking up group value        
        group_checked = True                                                   # if stays False then no group is used, should return none
 
        try:
            total_Hf += group_Hf[key]
        except:
            if return_missing_groups is False:                                 # return None if missing group value
                return None
            else:
                missing_groups_flag = True
                missing_groups.append(key)
    
    # apply ring correction
    try:
        ring_correction, missing_RC_flag = get_ring_correction_pathermo(m)
        if return_regardless_of_missing_RCs is False and missing_RC_flag is True:
            return None
    except:    # in case of error
        print("Error with ring correction function.")
        ring_correction = 0

    # Return total Hf
    if group_checked is True:
        if missing_groups_flag is False:
            return total_Hf + ring_correction 
        else:
            return missing_groups
        
    return None
