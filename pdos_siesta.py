import xml.etree.ElementTree as ET
import plotly.express as px
from statistics import mean

def read_siesta_pdos(project):
	tree = ET.parse(project)
	root = tree.getroot()

	for fermi in root.iter("fermi_energy"):
		fermi_energy_str = fermi.text.split()

	print('fermi: ' + fermi_energy_str[0])
	fermi_energy = float(fermi_energy_str[0])

	for energy in root.iter("energy_values"):
		energy_arr = energy.text.split()

	energy_arr = [float(i) - fermi_energy for i in energy_arr]
	energy_arr = [float(i) for i in energy_arr]

	occupations = []
	for orbital in root.iter('orbital'):
		dat = orbital.find('data')
		occupation = dat.text.split()
		occupation_float = [float(i) for i in occupation]
		occupations.append(occupation_float)

	return energy_arr, occupations

def read_siesta_pdos_twoSpecies(project, species1_upper_orbital_num):
	tree = ET.parse(project)
	root = tree.getroot()

	for fermi in root.iter("fermi_energy"):
		fermi_energy_str = fermi.text.split()

	print('fermi: ' + fermi_energy_str[0])
	fermi_energy = float(fermi_energy_str[0])

	for energy in root.iter("energy_values"):
		energy_arr = energy.text.split()

	energy_arr = [float(i) -fermi_energy for i in energy_arr]

	energy_arr = [float(i) for i in energy_arr]

	species1_occupations = []
	species2_occupations = []
	for orbital in root.iter('orbital'):
		index = orbital.get('index').split()
		dat = orbital.find('data')
		occupation = dat.text.split()
		occupation_float = [float(i) for i in occupation]
		if(int(index[0]) <= species1_upper_orbital_num):
			species1_occupations.append(occupation_float)
		else:
			species2_occupations.append(occupation_float)

	return energy_arr, species1_occupations, species2_occupations

def single_orbital_pdos(project, orbital_num):
	tree = ET.parse(project)
	root = tree.getroot()

	for energy in root.iter("energy_values"):
		energy_arr = energy.text.split()

	for fermi in root.iter("fermi_energy"):
		fermi_energy_str = fermi.text.split()

	fermi_energy = float(fermi_energy_str[0])
	energy_arr = [float(i) - fermi_energy for i in energy_arr]

	energy_arr = [float(i) for i in energy_arr]
	occupations = []
	for orbital in root.iter('orbital'):
		index = orbital.get('index').split()
		#print(index[0])
		if(index[0]==str(orbital_num)):
			n = orbital.get('n').split()
			l = orbital.get('l').split()
			m = orbital.get('m').split()
			p = orbital.get('P').split()
			print("n, l, m, polarized: "  + n[0] + " " + l[0] + " " + m[0] + " " + p[0])
			dat = orbital.find('data')
			occupation = dat.text.split()
			occupation_float = [float(i) for i in occupation]
			occupations.append(occupation_float)
			break

	return energy_arr, occupations

def orbital_by_species_pdos(project, species_name):
	tree = ET.parse(project)
	root = tree.getroot()

	for energy in root.iter("energy_values"):
		energy_arr = energy.text.split()

	for fermi in root.iter("fermi_energy"):
		fermi_energy_str = fermi.text.split()

	fermi_energy = float(fermi_energy_str[0])
	print('Fermi energy: ' + str(fermi_energy))
	energy_arr = [float(i) - fermi_energy for i in energy_arr]

	energy_arr = [float(i) for i in energy_arr]
	occupations = []
	for orbital in root.iter('orbital'):
		index = orbital.get('index').split()
		#print(index[0])
		species = orbital.get('species').split()
		#print(species)
		if(species[0]==str(species_name)):
			#print('Given Species Found')
			#n = orbital.get('n').split()
			#l = orbital.get('l').split()
			#m = orbital.get('m').split()
			#p = orbital.get('P').split()
			#print("n, l, m, polarized: "  + n[0] + " " + l[0] + " " + m[0] + " " + p[0])
			dat = orbital.find('data')
			occupation = dat.text.split()
			occupation_float = [float(i) for i in occupation]
			occupations.append(occupation_float)


	print('Total orbitals for given species (' + species_name +') = ' + str(len(occupations)))
	total_occupations = [0 for x in range(len(occupations[0]))]
	for i in range(len(occupations)):
		for j in range(len(occupations[i])):
			total_occupations[j] +=occupations[i][j]

	return energy_arr, total_occupations

def orbital_by_atom_index_pdos(project, atom_index):
	tree = ET.parse(project)
	root = tree.getroot()

	for energy in root.iter("energy_values"):
		energy_arr = energy.text.split()

	for fermi in root.iter("fermi_energy"):
		fermi_energy_str = fermi.text.split()

	fermi_energy = float(fermi_energy_str[0])
	print('Fermi energy: ' + str(fermi_energy))
	energy_arr = [float(i) - fermi_energy for i in energy_arr]

	energy_arr = [float(i) for i in energy_arr]
	occupations = []
	for orbital in root.iter('orbital'):
		index = orbital.get('atom_index').split()
		if(index[0]==str(atom_index)):
			#print('Given Species Found')
			dat = orbital.find('data')
			occupation = dat.text.split()
			occupation_float = [float(i) for i in occupation]
			occupations.append(occupation_float)


	print('Total orbitals for given index (' + str(atom_index) +') = ' + str(len(occupations)))
	total_occupations = [0 for x in range(len(occupations[0]))]
	for i in range(len(occupations)):
		for j in range(len(occupations[i])):
			total_occupations[j] +=occupations[i][j]

	return energy_arr, total_occupations	

def pdos_plot(project):
	pdos_values = read_siesta_pdos(project)
	energies = pdos_values[0]
	occupations = pdos_values[1]
	occupations_sum = [ sum(row[i] for row in occupations) for i in range(len(occupations[0])) ]

	fig = px.line(x=energies, y=occupations_sum)
	fig.update_layout(
    	title="DOS summation",
    	xaxis_title="Energy [eV]",
    	yaxis_title="PDOS",
    )
	fig.show()

def pdos_plot_twoSpecies(project, species1_upper_orbital_num, species1_name, species2_name):
	pdos_values = read_siesta_pdos_twoSpecies(project, species1_upper_orbital_num)
	energies = pdos_values[0]
	species1_occupations = pdos_values[1]
	species2_occupations = pdos_values[2]
	print(len(species2_occupations))
	species1_occupations_sum = [ sum(row[i] for row in species1_occupations) for i in range(len(species1_occupations[0])) ]
	species2_occupations_sum = [ sum(row[i] for row in species2_occupations) for i in range(len(species2_occupations[0])) ]
	print(len(species2_occupations_sum))

	fig = px.line(x=energies, y=species1_occupations_sum, labels={'y':species1_name})
	fig.add_scatter(x=energies, y=species2_occupations_sum, name=species2_name, mode='lines')
	fig.update_layout(
    	title="GroundState DOS Summation",
    	xaxis_title="Energy [eV]",
    	yaxis_title="PDOS",
    )
	fig.show()

def orbital_pdos_plot(project, orbital_num):
	pdos_values = single_orbital_pdos(project, orbital_num)
	energies = pdos_values[0]
	occupations = pdos_values[1]
	occupations_sum = [ sum(row[i] for row in occupations) for i in range(len(occupations[0])) ]

	fig = px.line(x=energies, y=occupations_sum)
	fig.update_layout(
    	title="DOS summation for orbital #" + str(orbital_num),
    	xaxis_title="Energy [eV]",
    	yaxis_title="PDOS",
    )
	fig.show()

def orbital_by_species_plot(project, species_name):
	pdos_values = orbital_by_species_pdos(project, species_name)
	energies = pdos_values[0]
	occupations = pdos_values[1]

	fig = px.line(x=energies, y=occupations)
	fig.update_layout(
    	title="DOS summation for:" + str(species_name),
    	xaxis_title="Energy [eV]",
    	yaxis_title="PDOS",
    )
	fig.show()

#returns the index of given array which contains value closest to given point
def get_closest(point, array):
	distance_from_array = [abs(i-point) for i in array]
	min_value = min(distance_from_array)
	return distance_from_array.index(min_value)

#find the s, p, or d center; say 'all' for entire sys, otherwise one atom_index given as int
def orbital_center(project, orb, atom_index):
	whole_system = False
	if str(atom_index)=='all':
		whole_system = True

	if orb =='s':
		l_value = 0
	elif orb =='p':
		l_value = 1
	elif orb =='d':
		l_value = 2
	tree = ET.parse(project)
	root = tree.getroot()

	for energy in root.iter("energy_values"):
		energy_arr = energy.text.split()

	for fermi in root.iter("fermi_energy"):
		fermi_energy_str = fermi.text.split()

	fermi_energy = float(fermi_energy_str[0])
	print('Fermi energy: ' + str(fermi_energy))
	energy_arr = [float(i) - fermi_energy for i in energy_arr]

	energy_arr = [float(i) for i in energy_arr]
	occupations = []
	max_occupation_energies = []
	mean_occupation_energies = []
	for orbital in root.iter('orbital'):
		l = orbital.get('l').split()
		index = orbital.get('atom_index').split()
		if whole_system:
			if(l[0]==str(l_value)):
				#print('Given Species Found')
				dat = orbital.find('data')
				occupation = dat.text.split()
				occupation_float = [float(i) for i in occupation]
				occupations.append(occupation_float)

				max_occupation = max(occupation_float)
				max_occupation_index = occupation_float.index(max_occupation)
				max_occupation_energies.append(energy_arr[max_occupation_index])

				mean_occupation_float = mean(occupation_float)
				mean_occupation_float_index = get_closest(mean_occupation_float, occupation_float)
				mean_occupation_energies.append(energy_arr[mean_occupation_float_index])
		else:
			if(l[0]==str(l_value) and index[0]==str(atom_index)):
				#print('Given Species Found')
				dat = orbital.find('data')
				occupation = dat.text.split()
				occupation_float = [float(i) for i in occupation]
				occupations.append(occupation_float)

				max_occupation = max(occupation_float)
				max_occupation_index = occupation_float.index(max_occupation)
				max_occupation_energies.append(energy_arr[max_occupation_index])

				mean_occupation_float = mean(occupation_float)
				mean_occupation_float_index = get_closest(mean_occupation_float, occupation_float)
				mean_occupation_energies.append(energy_arr[mean_occupation_float_index])


	#max_max_occupation_energy = max(max_occupation_energies)
	#mean_occupation_energy = mean(mean_occupation_energies)
	#print('Mode energy for ' +orb +' orbitals = ' + str(max_max_occupation_energy))
	#print('Mean energy for ' +orb +' orbitals = ' + str(mean_occupation_energy))

	total_occupations = [0 for x in range(len(occupations[0]))]
	for i in range(len(occupations)):
		for j in range(len(occupations[i])):
			total_occupations[j] +=occupations[i][j]

	max_occupation = max(total_occupations)
	max_occupation_index = total_occupations.index(max_occupation)
	max_occupation_energy = energy_arr[max_occupation_index]
	print('Energy of max occupation for ' +orb +' orbitals = ' + str(max_occupation_energy))
	return energy_arr, total_occupations
#make_pdos_plot("../../h2o.PDOS")
