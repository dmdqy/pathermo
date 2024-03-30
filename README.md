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
   >> Hf("CCC")
   0
```
If returns ```None```, it means the primary_groups file is missing one or more groups for your input molecule.


Optional flag: to see what groups are missing for an unsupported molecule, add the ```return_missing_groups = True ``` flag. It will return a list of missing groups if the molecule is unsupported.
```
   >> Hf("C#CO", return_missing_groups = True)
   ['CT 1 O', 'O 2 CT H']
```

In rare cases where a molecules is entirely made of groups that are not independent, it will return ```None``` and the ```return_missing_groups = True ``` flag won't take effect. In such special cases the data of the molecule can be added to the user_properties file if needed.
```
   >> Hf("N#CC#N", return_missing_groups = True)
   None
```


## References
Benson, Sidney William. "Thermochemical kinetics: methods for the estimation of thermochemical data and rate parameters." (No Title) (1976).
