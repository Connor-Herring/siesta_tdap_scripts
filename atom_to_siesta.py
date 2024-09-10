from ase import Atoms
from ase.io import read 


#takes input.fdf as a string, atoms object, number of species
def write_input_file(fdf_string, a, num_atoms):
	f = open("input.fdf", "w")
	f.write("SystemName " + str(a.symbols) + "\nSystemLabel siesta")
	f.write("\nNumberOfSpecies " + str(num_atoms) + "\nNumberOfAtoms " + str(len(a)))
	f.write(fdf_string)
	f.close()

#finds number identifier of a repeated atom (used in atoms_to_siesta())
def find_location(number, arr):
	for x in range(0, len(arr)):
		if(arr[x] == number):
			return (x+1)

#main script
def atoms_to_siesta(a):
	f = open("coordinates.data", "w") #create struct_in file
	f2 = open("label.data", "w") #chemical species label file
	f3 = open("lattice.data", "w") #holds lattice vector
	numAtoms = len(a) #number of atoms in atoms object
	print(numAtoms)
	atomicNum = a.get_atomic_numbers() #array of atomic numbers
	count = 0 #keeps track of how many unique atoms there are
	temp = [0 for x in range(numAtoms)] #used to check if an atom is unique
	atomCount = [0 for x in range(numAtoms)] #each list element corresponds to unique atom number
	atomPos = a.get_scaled_positions() #holds position of each atom 
	label = a.get_chemical_symbols() #holds chemical symbols
	unique_number_list = []

	#write the cell to lattice file-----------------
	for i in range(0, 3):
		for j in range(0,3):
			f3.write(str(a.cell[i][j]) + " ")
		f3.write("\n")
	f3.write("\n")
	#-------------------------------------------

	#loop through each atom in Atoms
	for x in range(0, numAtoms):
		for y in range(0,3):
			f.write(str(round(atomPos[x][y],6)) + " ")

		if (atomicNum[x] in temp): #if atom is repeated (not unique), find its indentifier number
			location = find_location(atomicNum[x], unique_number_list)
			atomCount[x] = location
		else:
			unique_number_list.append(atomicNum[x])
			count+=1 #increment unique atom count
			f2.write(str(count) + " " + str(atomicNum[x]) + " " + str(label[x]) + "\n") #write to species label file
			atomCount[x] = count

		temp[x] = atomicNum[x]
		f.write("  " + str(atomCount[x]) +  " #" + str(x+1))
		f.write("\n")
	f. close()
	return count




