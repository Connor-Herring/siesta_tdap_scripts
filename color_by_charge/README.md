These scripts are used after running a TDAP calculation to color each atom by its charge. Make sure your calculation wrote Hirshfeld populations before proceeding. 

### Step 1. Generate Mol2 Files
##### In a notebook run something like the following to create a mol2 file for every nth geometry. I would recommend creating at least 60 so the final movie has a reasonable number of frames. This example code takes steps of 50 with a total of 3000 MD steps in the simulation. 

```
from siesta_to_mol2 import *
for i in range (1, 3001, 50):
    data = get_coords_and_charge(projectDirectoryName, i)
    write_mol2(data[0], data[1], str(i) + '.mol2', '../../../calculations/Cu_54_Ru_N2_vertical_2.5eV_0.1amplitude_fixedPositions/mol2s/')
```

### Step 2. Render
##### Open VMD and render the geometries using the following two commands (you will need to adjust for your paths). The second path should be the path to the directory containing the mol2s created in step 1. 

```
source C:/Users/Dell/Desktop/Lab/vmdScripts/charge_test.tcl
main_loop C:/Users/Dell/Desktop/Lab/DFT/Calculations/Cu_54_Re_N2_vertical_2.5eV_0.1amplitude_fixedPositions/mol2s/
```

### Step 3. Make Movie
##### Stitch together rendered images to a movie. 

```
python make_movie.py
```
