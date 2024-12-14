import shutil

calculation_path = '../../../calculations/'

#returns a list of atomic coordinates and corresponding charges
def get_coords_and_charge(project, md_step):
	file =open(calculation_path + project + '/out', 'r')
	Lines = file.readlines()
	current_line = 0
	counter_check =0
	bulk_charge_at_time = []
	time_step = 0.0242
	current_md_step = 0

	#only grabs the last hirshfeld atomic charge before next MD step, hence convoluted nested loops :)
	for k in range(0, len(Lines)):
		cut = Lines[k].split()
		current_line = k
		if(len(cut)>1):
			if cut[0]=='NumberOfAtoms':
				num_atoms = int(cut[1])
				print('number of atoms found: ' + str(num_atoms))
			if(cut[0] == 'Begin' and cut[1] =='MD'):
				if cut[4]==str(md_step):
					print(cut)
					while(True):
						try:
							new_cut = Lines[current_line].split()
							if(len(new_cut) >1 and new_cut[0] =='Hirshfeld' and new_cut[1]=='Net'):
								charges = Lines[current_line+2:current_line+num_atoms+2]
								charges_to_return = []
								for elem in charges:
									split_elem = elem.split()
									charges_to_return.append(split_elem)

								current_line += num_atoms

							elif(len(new_cut) >1 and new_cut[0]=='*' and new_cut[1]=='Maximum'):
								break
							current_line +=1
							k=current_line
						except Exception as e:
							print(e)
							print(new_cut)
							break


				current_md_step +=1

			if cut[0]=='siesta:' and cut[1]=='Atomic' and cut[2]=='coor.' and current_md_step==md_step:
				coords = Lines[current_line+1:current_line + num_atoms+1]
				coords_to_return = []
				for elem in coords:
					split_elem = elem.split()
					coords_to_return.append(split_elem[1:4])

	return coords_to_return, charges_to_return
	

#writes coordinates and charges to mol2 file 
def write_mol2(coords, charges, filename, path_to_move):
	num_atoms = len(coords)
	atoms_dict = {}

	with open(filename, 'w') as file:
		file.write('@<TRIPOS>MOLECULE\n' + 'test\n'+ '   ' + str(num_atoms) + ' ' + str(num_atoms) + ' 0 0 0\n' + 'SMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n')

		for i in range(num_atoms):
			if charges[i][2] in atoms_dict.keys():
				curr_count = atoms_dict[charges[i][2]]
				atoms_dict.update({str(charges[i][2]): curr_count + 1})#add 1 to count
				atom_id = curr_count+1
			else:
				atoms_dict.update({str(charges[i][2]): 1})#create atom and its count in dictionary 
				atom_id = 1

			file.write('      ' +str(i+1) +' ' + str(charges[i][2]) +str(atom_id) +' ' + str(coords[i][0]) +' ' + str(coords[i][1]) +' ' + str(coords[i][2]) + ' ' +str(charges[i][2]) + ' 1 test ' + str(charges[i][1])+'\n')

		file.write('@<TRIPOS>BOND\n')
		for i in range(num_atoms):
			file.write('1    1    2 ar\n')

	shutil.move(filename, path_to_move +filename)

