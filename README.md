# pathermo
This code is designed for predicting the standard enthalpy of formation using the Benson group-increment theory. 'pathermo' is a sanitized version of 'pgthermo', originally developed and utilized internally within the Broadbelt group at Northwestern University. 'pathermo' includes a restricted selection of published Benson group values. Users have the flexibility to incorporate their own group values for their specific requirements.




## Prerequisites 

rdkit >= 2023.03.2

pandas



## Installation 
1. Navigate to the installation folder
   ```
   cd ~/GitHub-Repos
   ```
2. Clone the repo, navigate to the folder
   ```
   git clone https://github.com/dmdqy/pathermo.git
   cd ~/GitHub-Repos/pathermo
   ```
3. Install using pip
   ```
   pip install -e .
   ```

## Update Program
In case the code is updated, you may use git pull to update yours to the latest version.

1. Navigate to the pgthermo folder
   ```   
   cd ~/GitHub-Repos/pathermo
   ```
2. Pull. If you used the -e flag as recommended when installing, it will update automatically after pull.
   ```
   git pull
   ```



## Usage
pgthermo takes SMILES of a molecule or radical as input, and returns the gas-phase data of standard enthalpy of formation.
```
   >> from pathermo.properties import Hf
   >> Hf("C1CCCCC1")
   -29.58
```
If returns ```None```, it means the primary_groups file is missing one or more groups for your input molecule.


Optional flag: to see what groups are missing for an unsupported molecule, add the ```return_missing_groups = True ``` flag. It will return a list of missing groups if the molecule is unsupported.
```
   >> Hf("C#CO", return_missing_groups = True)
   ['CT 1 H', 'CT 1 O', 'O 2 CT H']
```

In rare cases where a molecules is entirely made of groups that are not independent, it will return ```None``` and the ```return_missing_groups = True ``` flag won't take effect. In such special cases the data of the molecule can be added to the user_properties file if needed.
```
   >> Hf("N#CC#N", return_missing_groups = True)
   None
```


## References
[1] Benson, Sidney William. "Thermochemical kinetics: methods for the estimation of thermochemical data and rate parameters." (1976).  
[2] Eigenmann, H. K., D. M. Golden, and S. W. Benson. "Revised group additivity parameters for the enthalpies of formation of oxygen-containing organic compounds." The Journal of Physical Chemistry 77, no. 13 (1973): 1687-1691.  
[3] Patai, Saul, Zvi Rappoport, and Charles James Matthew Stirling. "The chemistry of sulphones and sulphoxides." (1988).  
[4] Gillis, Ryan J., and William H. Green. "Thermochemistry prediction and automatic reaction mechanism generation for oxygenated sulfur systems: A case study of dimethyl sulfide oxidation." ChemSystemsChem 2, no. 4 (2020): e1900051.  
[5] Holmes, John L., and Christiane Aubry. "Group additivity values for estimating the enthalpy of formation of organic compounds: an update and reappraisal. 2. C, H, N, O, S, and halogens." The Journal of Physical Chemistry A 116, no. 26 (2012): 7196-7209.  
[6] Ashcraft, Robert W., and William H. Green. "Thermochemical properties and group values for nitrogen-containing molecules." The Journal of Physical Chemistry A 112, no. 38 (2008): 9144-9152.  
[7] Domalski, Eugene S., and Elizabeth D. Hearing. "Estimation of the thermodynamic properties of C‐H‐N‐O‐S‐halogen compounds at 298.15 K." Journal of Physical and Chemical Reference Data 22, no. 4 (1993): 805-1159.  
[8] P.J. Linstrom and W.G. Mallard, Eds., NIST Chemistry WebBook, NIST Standard Reference Database Number 69, National Institute of Standards and Technology, Gaithersburg MD, 20899, https://doi.org/10.18434/T4D303, (retrieved March 29, 2024).  
