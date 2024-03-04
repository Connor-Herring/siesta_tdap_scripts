import numpy as np
import numpy.fft as npf 
import matplotlib.pyplot as plt
import math
from statistics import mean
import plotly.express as px
import plotly.graph_objects as go
from hirshfeld_charge import *

dt = 0.0242 #in fs (same as 0.5 hbar/Ry)
#dt = 1 #in fs
#dt = 0.5 #in hbar/ry
pi = 3.14159
hbarev =  6.58211928000000e-16

#a of form [[x,y,z], [x,y,z],.....] i.e. components of dipole for each time step
def plot_fft(a):
	x_dipole = [] #hold x component of dipole
	y_dipole = [] #hold y component of dipole
	z_dipole = [] #hold z component of dipole
	x_fft = [] #fft of x component
	y_fft = [] #fft of y component
	z_fft = [] #fft of z component
	imag_sum = [] #sum of imaginary parts of x,y,z fft

	#throw components into corresponding arrays
	for i in range(len(a)):
		x_dipole.append(a[i][0])
		y_dipole.append(a[i][1])
		z_dipole.append(a[i][2])

	#take fft of 3 comps
	x_fft = npf.fft(x_dipole, n=None, axis=-1, norm=None)
	y_fft = npf.fft(y_dipole, n=None, axis=-1, norm=None)
	z_fft = npf.fft(z_dipole, n=None, axis=-1, norm=None)

	for j in range(0, len(a)):
		imag_sum.append(x_fft.imag[j] + y_fft.imag[j] + z_fft.imag[j])

	freq = npf.fftfreq(len(a), d=0.0242) #units of 10^15 Hz
	freq_converted = [x/.2417989 for x in freq] #converet to ev
	#dw= 2.0*pi*hbarev/(dt*1E-15)/(len(imag_sum)*2)
	#x=np.linspace(1.*dw,(len(imag_sum)-1)*dw,(len(imag_sum)-1))
	#plt.plot(x, imag_sum[1:])
	plt.plot(freq_converted, imag_sum)
	#plt.xticks(np.arange(min(freq), max(freq)+1, 0.1))
	plt.xlim(0,10)
	plt.ylim(0,700)
	#plt.ylim(0,.003)
	plt.xlabel('fft [ev]')
	plt.show()

def plot_dipole(proj1, proj2, proj3):
	f1 = open(calculation_path + proj1 + '/dipole.vs.time', 'r')
	f2 = open(calculation_path + proj2 + '/dipole.vs.time', 'r')
	f3 = open(calculation_path + proj3 + '/dipole.vs.time', 'r')
	Lines1 = f1.readlines()
	Lines2 = f2.readlines()
	Lines3 = f3.readlines() 
	Lines_list = [Lines1, Lines2, Lines3]
	data_points = len(Lines_list[0]) #start by assuing first file is shortest

	#find the shortest file incase the lengths don't match to avoid issues
	for i in range(1, 3):
		if(len(Lines_list[i])< data_points):
			data_points = len(Lines_list[i])

	print(data_points)
	time = []
	dipole_magnitude = []
	dipole_components = []

	for i in range(0,3):
		#Lines = files[i].readlines()
		array_loc = 0
		for line in Lines_list[i]:
			if(array_loc == data_points): #to ensure equal data points in all dipole lists
				break
			cut = line.split()

			if(i==0):
				t = float(cut[0])*.0242  #.0242 is time in fs for each step
				time.append(t)
				dipole_components.append([]) #build list once [[x,y,z], [x,y,z]....etc]

			value = float(cut[i+1])
			dipole_components[array_loc].append(value) #append all the x's, then all y's, then all z's
			array_loc +=1

	#calculate magnitude of dipole
	for i in range(0, len(dipole_components)):
		x = dipole_components[i][0]
		y = dipole_components[i][1]
		z = dipole_components[i][2]
		mag = math.sqrt((x*x) + (y*y) + (z*z))
		dipole_magnitude.append(mag)

	plt.xlabel('t (fs)')
	plt.ylabel('dipole ')
	plt.plot(time, dipole_magnitude)
	plt.show()
	#plot_fft(dipole_components)

def get_dipole(path):
	f = open(path + '/dipole.vs.time', 'r')
	Lines = f.readlines()

	time = []
	dipole_magnitude = []

	for line in Lines:
		cut = line.split()
		time.append(float(cut[0])*.0242)
		x = float(cut[1])
		y = float(cut[2])
		z = float(cut[3])
		dipole_magnitude.append(math.sqrt(x**2 + y**2 + z**2))


	return [time, dipole_magnitude]

def get_energy(path):
	f = open(path + '/energy.vs.time', 'r')

	Lines = f.readlines()
	time = []
	total_energy = []

	for line in Lines:
		cut = line.split()
		if(cut[0] != '#'):
			total_energy.append(float(cut[1]))
			t = float(cut[0])*.0242  #.0242 is time in fs for each step
			time.append(t)

	return [time, total_energy]

def plot_energy(path):
	#f = open('../../calculations/' + proj1 + '/siesta.MDE', 'r')
	data = get_energy(path)
	time = data[0]
	total_energy = data[1]

	fig = px.line(x=time, y=total_energy)
	fig.update_layout(
	    	xaxis_title="Time (fs)",
	    	yaxis_title="E (ev)",
	    )
	fig.show()
	print('initial energy: ' + str(total_energy[0]))
	print('final energy: ' + str(total_energy[-1]))

	
#used for ANI file ****atom1 number must precede atom2**** result is in Ang
def get_bondlength(path, atom1, atom2):
	time_step = dt
	filename = path + '/siesta.ANI'
	f = open(filename, 'r')
	Lines = f.readlines()	
	current_line = 0
	bond_length = []
	num_atoms = int(Lines[0])
	atom1_symbol = str(Lines[atom1 + 1].split()[0])
	atom2_symbol = str(Lines[atom2 + 1].split()[0])
	print(atom1_symbol + " " + atom2_symbol)
	print(Lines[atom1 + 1].split())
	print(Lines[atom2 + 1].split())
	num_diff = atom2 - atom1
	current_line += atom1

	for i in range(len(Lines)):
		if current_line>=len(Lines):
			break
		cut = Lines[current_line].split()
		if(len(cut)==1):
			if(int(cut[0])==num_atoms):
				atom1_pos = Lines[current_line + atom1 +1].split()
				atom2_pos = Lines[current_line +atom2 +1].split()
				a = float(atom1_pos[1]) - float(atom2_pos[1])
				b = float(atom1_pos[2]) - float(atom2_pos[2])
				c = float(atom1_pos[3]) - float(atom2_pos[3])
				length = math.sqrt((a*a)  +(b*b) + (c*c))
				bond_length.append(length)
		current_line +=1

	steps = list(range(0, len(bond_length)))
	print('steps: ' + str(len(steps)))
	fs_steps = [time_step*i for i in steps]

	return [fs_steps, bond_length]

def plot_bondlength(path, atom1, atom2):
	data = get_bondlength(path, atom1, atom2)
	fs_steps = data[0]
	bond_length = data[1]

	fig = px.line(x=fs_steps, y=bond_length)
	fig.update_layout(
			title="Bond length vs Time",
	    	xaxis_title="Time (fs)",
	    	yaxis_title="Bond length (Ang)",
	    )
	fig.show()
	initial_avg_bondlength = sum(bond_length[0:10])/10
	print('initial bond length: ' + str(initial_avg_bondlength))
	final_avg_bondlength = sum(bond_length[-10:])/10
	print("final avg bond length: " + str(final_avg_bondlength))
	delta_bond = final_avg_bondlength - initial_avg_bondlength
	print('delta: ' + str(delta_bond))

#Returns time point (in steps not fs) in simulation where atom1 and atom2 bond length is >= given one 
def get_time_at_bondlength(path, atom1, atom2, bondlength_to_find):
	time_step = 0.0242
	filename = path + '/siesta.ANI'
	f = open(filename, 'r')
	Lines = f.readlines()	
	current_line = 0
	time_step_count = 0
	initial_bondlength = 0
	num_atoms = int(Lines[0])
	atom1_symbol = str(Lines[atom1 + 1].split()[0])
	atom2_symbol = str(Lines[atom2 + 1].split()[0])
	print(atom1_symbol + " " + atom2_symbol)
	print(Lines[atom1 + 1].split())
	print(Lines[atom2 + 1].split())
	num_diff = atom2 - atom1
	current_line += atom1

	for i in range(len(Lines)):
		if current_line>=len(Lines):
			break
		cut = Lines[current_line].split()
		if(len(cut)==1):
			if(int(cut[0])==num_atoms):
				atom1_pos = Lines[current_line + atom1 +1].split()
				atom2_pos = Lines[current_line +atom2 +1].split()
				a = float(atom1_pos[1]) - float(atom2_pos[1])
				b = float(atom1_pos[2]) - float(atom2_pos[2])
				c = float(atom1_pos[3]) - float(atom2_pos[3])
				current_length = math.sqrt((a*a)  +(b*b) + (c*c))

				if time_step_count ==0:
					initial_bondlength = current_length

				bond_length_delta = current_length - initial_bondlength
				time_step_count +=1
				time_fs = time_step_count*.0242
				
				if current_length >= bondlength_to_find:
					print('initial_bond length: ' + str(initial_bondlength))
					print('current bond length: ' + str(current_length))
					print(str(time_step_count) + ' or ' + str(time_step_count*.0242) + ' fs')
					break

				#if bond_length_delta >= bondlength_to_find:
					#print('initial_bond length: ' + str(initial_bondlength))
					#print('current bond length: ' + str(current_length))
					#print(str(time_step_count) + ' or ' + str(time_step_count*.0242) + ' fs')
					#break
		current_line +=1

	return [time_fs, time_step_count]

def get_bondlength_delta(path, atom1, atom2):
	data = get_bondlength(path, atom1, atom2)
	bl_delta = [i - data[1][0] for i in data[1]]

	return [data[0], bl_delta]

def get_avg_time_at_bondlength(path, main_atom, atom_start_number, atom_end_number, bondlength_to_find):

	fs_times = []
	time_steps = []
	for atom in range(atom_start_number, atom_end_number+1):
		data = get_time_at_bondlength(path, main_atom, atom, bondlength_to_find)
		fs_times.append(data[0])
		time_steps.append(data[1])

	avg_fs_time = mean(fs_times)
	avg_time_step = int(mean(time_steps))

	print(fs_times)
	print('----------------avg trans state time: ' + str(avg_fs_time))
	return [avg_fs_time, avg_time_step]

#plots bondlength in Ang between atom 1 and 2 from siesta out file
def plot_bondlength_from_out(path, atom1, atom2):
	time_step = 0.0242
	filename = path + '/out'
	f = open(filename, 'r')
	Lines = f.readlines()	
	current_line = 0
	bond_length = []
	num_atoms = 68

	for i in range(len(Lines)):
		if current_line>=len(Lines):
			break
		cut = Lines[current_line].split()
		if(len(cut)==4 and cut[-1]=='(Ang):'):
			atom1_pos = Lines[current_line + atom1].split()
			atom2_pos = Lines[current_line + atom2].split()
			a = float(atom1_pos[1]) - float(atom2_pos[1])
			b = float(atom1_pos[2]) - float(atom2_pos[2])
			c = float(atom1_pos[3]) - float(atom2_pos[3])
			length = math.sqrt((a*a)  +(b*b) + (c*c))
			bond_length.append(length)
			current_line +=num_atoms
		current_line +=1

	steps = list(range(0, len(bond_length)))
	fs_steps = [time_step*i for i in steps]
	fig = px.line(x=fs_steps, y=bond_length)
	fig.update_layout(
			title="bond length",
	    	xaxis_title="Time (fs)",
	    	yaxis_title="BL (Ang)",
	    )
	fig.show()

def summary(path):
	filename = path + '/out'
	f = open(filename, 'r')
	Lines = f.readlines()	
	current_line = 0
	md_step = 0
	max_scf_steps = 0
	current_line = 0
	md_steps = []
	scf_count_per_step = []

	for i in range(len(Lines)):
		if current_line>=len(Lines):
					break

		line = Lines[current_line]
		if "Dump of input data file" in line: #if we are in input dump section iterate through it and print things of interest
			print("-----Important input data-------")
			while True:
				line = Lines[current_line]
				if "NumberOfSpecies" in line:
					num_species = line.split()[1]
				elif "NumberOfAtoms" in line:
					num_atoms = line.split()[1]
				elif "SpinPolarized" in line:
					if(len(line.split()) <3):
						spin_polarized = line.split()[1]
				elif "TD.ntime" in line:
					ntime = line.split()[1]
				elif "MaxSCFIterations" in line:
					max_scf_steps = line.split()[1]
				elif "End of input data file" in line:
					break
				current_line +=1

		
		if "Begin MD step" in line:
			#print(line.split())
			md_steps.append(line.split()[4])
			while True:
				try:
					line = Lines[current_line+1]
				except:
					break
				if "siesta" in line and len(line.split()) ==7:
					if line.split()[1] != 'iscf':
						scf_step = line.split()[1]
						#print(line)
				elif "Begin MD step" in line:
					current_line -=2
					scf_count_per_step.append(scf_step)
					break
				current_line +=1
	
		current_line +=1

	print("NumberOfSpecies: " + num_species)
	print("NumberOfAtoms: " + num_species)
	print("SpinPolarized: " + spin_polarized)
	print("TD.ntime: " + ntime)
	print("max_scf_steps: " + max_scf_steps)
	print("--------------------")
	if(len(md_steps) != len(scf_count_per_step)):
		md_steps.remove(md_steps[-1])
	fig = px.line(x=md_steps, y=scf_count_per_step)
	fig.update_layout(
	    	xaxis_title="md step",
	    	yaxis_title="scf steps",
	    )
	fig.show()

def check_converge_old(path):
	filename = path + '/out'
	f = open(filename, 'r')
	Lines = f.readlines()	
	current_line = 0
	dDmax = []
	Ef = []
	iscf =[]

	for i in range(len(Lines)):
		if current_line>=len(Lines):
			break
		cut = Lines[current_line].split()
		if(len(cut)>1):
			if(cut[0]=='siesta:' and cut[1]=='iscf'): #found scf loop for a given CG move, now loop to end
				while(current_line+1 <len(Lines)):
					current_line +=1
					scf_line = Lines[current_line]
					scf_cut = scf_line.split()
					if(len(scf_cut) <1):
						iscf.append(last_line[1]) #how many iscf iterations we hit
						break
					if(scf_cut[0] == 'siesta:'):
						try:
							dDmax.append(float(scf_cut[5]))
							Ef.append(float(scf_cut[6]))
						except:
							print('bad line: ' + str(scf_cut))
					last_line = scf_cut		
		current_line +=1
	print(dDmax[0])
	steps = list(range(0, len(dDmax)))
	fig = px.line(x=steps, y=dDmax)
	fig.update_layout(
			title="dDmax",
	    	xaxis_title="steps",
	    	yaxis_title="dDmax",
	    )
	fig.show()
	fig2 = px.line(x=steps, y=Ef)
	fig2.update_layout(
			title="Ef",
	    	xaxis_title="steps",
	    	yaxis_title="Ef",
	    )
	fig2.show()
	iscf_steps =list(range(0, len(iscf)))
	fig3 = px.line(x=iscf_steps, y=iscf)
	fig3.update_layout(
			title="iscf",
	    	xaxis_title="steps",
	    	yaxis_title="iscf",
	    )
	fig3.show()

#get the net spin (up - down) at each step 
def get_spin(path, atom_number):
	file =open(path + '/out', 'r')
	Lines = file.readlines()
	atom_spin  =[]
	time_step = 0.0242

	for k in range(0, len(Lines)):
		cut = Lines[k].split()
		current_line = k
		if(len(cut)>1):
			if(cut[0] == 'Begin' and cut[1] =='MD'):
				current_spin_set = False
				while(True):
					new_cut = Lines[current_line].split()
					if(len(new_cut) >2 and new_cut[0] =='Spin-projected' and new_cut[1] == 'Hirshfeld'):
						spin_string = str(Lines[current_line+atom_number+1]).split()
						#print(spin_string)
						if(len(atom_spin)<1):
							atom_letter = str(spin_string[3])
							print(atom_letter)

						current_spin = abs(float(spin_string[1]) - float(spin_string[2]))
						current_spin_set = True

					elif(current_spin_set and len(new_cut) >1 and new_cut[0]=='*' and new_cut[1]=='Maximum'):
						atom_spin.append(current_spin)
						break
					current_line +=1
					k +=1
	print('length of atom spin: ' + str(len(atom_spin)))
	atom_spin.remove(atom_spin[0])
	steps = list(range(0, len(atom_spin)))
	fs_time = [time_step*i for i in steps]
	return[fs_time, atom_spin]

#get the net spin (up - down) at each step from start atom to end atom
def get_spin_over_range(path, atom_start_number, atom_end_number):
	file =open(path + '/out', 'r')
	Lines = file.readlines()
	total_spin  =[]
	time_step = 0.0242

	for k in range(0, len(Lines)):
		cut = Lines[k].split()
		current_line = k
		if(len(cut)>1):
			if(cut[0] == 'Begin' and cut[1] =='MD'):
				current_spin_set = False
				while(True):
					new_cut = Lines[current_line].split()
					if(len(new_cut) >2 and new_cut[0] =='Spin-projected' and new_cut[1] == 'Hirshfeld'):
						sys_spin = 0
						for line_number in range(current_line + atom_start_number +1, atom_end_number):
							spin_string = str(Lines[line_number]).split()
							#print(spin_string)
							if(len(atom_spin)<1):
								atom_letter = str(spin_string[3])
								print(atom_letter)

							current_spin = abs(float(spin_string[1]) - float(spin_string[2]))
							current_spin_set = True
							sys_spin +=current_spin

					elif(current_spin_set and len(new_cut) >1 and new_cut[0]=='*' and new_cut[1]=='Maximum'):
						total_spin.append(sys_spin)
						break
					current_line +=1
					k +=1
	print('length of atom spin: ' + str(len(total_spin)))
	total_spin.remove(total_spin[0])
	steps = list(range(0, len(total_spin)))
	fs_time = [time_step*i for i in steps]
	return[fs_time, total_spin]

def plot_spin(path, atom_number):
	data = get_spin(path, atom_number)

	fig = px.line(x=data[0], y=data[1])
	fig.update_layout(
			title="Spin for atom: " + str(atom_number),
	    	xaxis_title="Time (fs)",
	    	yaxis_title="Net Spin",
	    )

	fig.show()

def plot_spin_over_range(path, atom_start_number, atom_end_number):
	data = get_spin_over_range(path, atom_start_number, atom_end_number)

	fig = px.line(x=data[0], y=data[1])
	fig.update_layout(
			title="Spin for atom: " + str(atom_start_number) + ' to ' + str(atom_end_number),
	    	xaxis_title="Time (fs)",
	    	yaxis_title="Net Spin",
	    )

	fig.show()

def autoNags_charge(path):
	time_step = 3000*.0242/60
	gs_charge = []
	na_charge = []
	na_filename = path + '/NA_charge'
	gs_filename = path + '/GS_charge'
	gs_file = open(gs_filename, 'r')
	gs_Lines = gs_file.readlines()
	na_file = open(na_filename, 'r')
	na_Lines = na_file.readlines()

	for i in range(len(gs_Lines)):
		gs_charge.append(float(gs_Lines[i]))
		na_charge.append(float(na_Lines[i]))

	steps = list(range(len(gs_charge)))
	fs_time = [time_step*i for i in steps]

	fig = go.Figure()
	fig.add_trace(go.Scatter(x=fs_time, y=gs_charge,
                    mode='lines',
                    name='groundstate'))
	fig.add_trace(go.Scatter(x=fs_time, y=na_charge,
                    mode='lines',
                    name='non adiabatic'))
	fig.update_layout(
			title="autoNags_charge",
	    	xaxis_title="Time (fs)",
	    	yaxis_title="Charge (e)",
	    )
	fig.show()

	charge_diff = [na_charge[i] - gs_charge[i] for i in range(len(na_charge))]
	fig2 =px.line(x=fs_time, y=charge_diff)
	fig2.update_layout(
			title="charge_diff",
	    	xaxis_title="Time (fs)",
	    	yaxis_title="charge (e)",
	    )
	fig2.show()


def get_atomic_coordinates(path, atom):
	filename = path + '/siesta.ANI'
	f = open(filename, 'r')
	Lines = f.readlines()	
	num_atoms = int(Lines[0])
	current_line = 0
	atom_coords = []

	for i in range(len(Lines)):
		cut = Lines[i].split()
		if len(cut) >0 and cut[0]==str(num_atoms):
			coords = Lines[i + atom +1].split()
			if i==0:
				print('Atom: ' + str(coords[0]))
			atom_coords.append(coords[1:4]) #only append the coords, not atom symbol

	atom_coords_float = [[float(x) for x in sub_list] for sub_list in atom_coords]#convert to floats
	return(atom_coords_float)

def distance_between(atom1_pos, atom2_pos):
	a = float(atom1_pos[1]) - float(atom2_pos[1])
	b = float(atom1_pos[2]) - float(atom2_pos[2])
	c = float(atom1_pos[3]) - float(atom2_pos[3])
	length = math.sqrt((a*a)  +(b*b) + (c*c))

	return length

#atom_list given as: [56, 57, 58] for example
def calculate_dipole(path, atoms_list):
	all_atoms_charge = []
	all_atoms_coordinates = []

	for atom in atoms_list:
		charge_data = get_charge_single_outFile(project, atom)
		all_atoms_charge.append(charge_data[1])

		coord_data = get_atomic_coordinates(project, atom)
		coord_data.remove(coord_data[0])
		all_atoms_coordinates.append(coord_data)

	dipole = []
	for i in range(len(all_atoms_coordinates[0])):
		component1 = (all_atoms_charge[0][i] + all_atoms_charge[1][i])*distance_between(all_atoms_coordinates[0][i], all_atoms_coordinates[1][i])
		component2 = (all_atoms_charge[0][i] + all_atoms_charge[2][i])*distance_between(all_atoms_coordinates[0][i], all_atoms_coordinates[2][i])
		component3 = (all_atoms_charge[1][i] + all_atoms_charge[2][i])*distance_between(all_atoms_coordinates[1][i], all_atoms_coordinates[2][i])
		dipole.append(component1 + component2 + component3)

	return dipole

#get the bond angle between 3 atoms, arguments [x,y,z], [x,y,z], [x,y,z]
def get_bond_angle(atom1_coords, center_atom_coords, atom3_coords):
	# Convert coordinates to numpy arrays
    atom1 = np.array(atom1_coords)
    atom2 = np.array(center_atom_coords)
    atom3 = np.array(atom3_coords)
    
    # Vectors between atoms
    vector1 = atom1 - atom2
    vector2 = atom3 - atom2
    
    # Calculate dot product and magnitudes
    dot_product = np.dot(vector1, vector2)
    magnitude1 = np.linalg.norm(vector1)
    magnitude2 = np.linalg.norm(vector2)
    
    # Calculate cosine of angle
    cosine_angle = dot_product / (magnitude1 * magnitude2)
    
    # Calculate angle in radians
    angle_radians = np.arccos(cosine_angle)
    
    # Convert angle to degrees
    angle_degrees = np.degrees(angle_radians)
    
    return angle_degrees
