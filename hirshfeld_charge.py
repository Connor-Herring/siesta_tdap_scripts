import numpy as np
import numpy.fft as npf 
import matplotlib.pyplot as plt
import math
from statistics import mean
import plotly.express as px
import plotly.graph_objects as go


#iterates from 0 to number of out files and passes each file onto get_charge which returns the charge of given atom. These values are stored on atom_charge[]
#atom_charge[] is then ploted 
def plot_charge(path, num_files, atom_number):
	count = 0
	fill_num = 0 #used for correct number of 0s in out filenames 
	atom_charge = []

	for x in range(0, num_files+1):
		if(count<10):
			filename = path + '/out.000' + str(fill_num) + str(count)
			f = open(filename, 'r')
			charge_data = get_charge(f, atom_number)
			atom_charge.append(float(charge_data[1]))
			print('plotting charge for: ' + charge_data[2])
		else: 
			filename = path + '/out.000' + str(count)
			f = open(filename, 'r')
			charge_data = get_charge(f, atom_number)
			try:
				atom_charge.append(float(charge_data[1]))
			except:
				print('error')
		count +=1

	steps = list(range(0, len(atom_charge)))
	#o_charge = [-.23]*len(atom_charge)
	#ob_charge = [0.142]*len(atom_charge)
	#oc_charge = [0.097]*len(atom_charge)
	plt.plot(steps, atom_charge)
	#plt.plot(steps, o_charge, 'g')
	#plt.plot(steps, ob_charge, 'r')
	#plt.plot(steps, oc_charge, 'm')
	plt.xlabel('steps')
	plt.ylabel('charge (e)')
	plt.show()
	print('intial charge: ' + str(atom_charge[0]))
	max_charge = avg_max_charge(atom_charge)
	print("rolling avg max charge: " + str(max_charge))
	min_charge = avg_min_charge(atom_charge)
	print("rolling avg min charge: " + str(min_charge))
	print(len(atom_charge[-10:]))
	avg = sum(atom_charge[-10:])/10
	print("final charge: " + str(avg))
	

#searches through out file to find the charge of a given atom, returns charge value and atom name
def get_charge(file, atom_number):
	Lines = file.readlines()
	current_line = 0

	for line in Lines:
		cut = line.split()
		if(len(cut)>1):
			if(cut[0] =='Hirshfeld'):
				charge_string = str(Lines[current_line+atom_number+1]).split()
				return charge_string
		current_line +=1

def avg_max_charge(atom_charge_array):
	max_charge = max(atom_charge_array)
	max_index = atom_charge_array.index(max_charge)

	rolling_avg = sum(atom_charge_array[max_index-10:max_index+10])/21
	return rolling_avg

def avg_min_charge(atom_charge_array):
	min_charge = min(atom_charge_array)
	min_index = atom_charge_array.index(min_charge)

	rolling_avg = sum(atom_charge_array[min_index-10:min_index+10])/21
	return rolling_avg

#-----------------------------------------------------------------------------------------------------------------------------------------#
#if you only have 1 out file with hirshfeld being written per geometry, use this
def get_charge_single_outFile(path, atom_number):
	file =open(path + '/out', 'r')
	Lines = file.readlines()
	current_line = 0
	counter_check =0
	atom_charge = []
	time_step = 0.0242
	current_charge = 0

	#only grabs the last hirshfeld atomic charge before next MD step, hence convoluted nested loops :)
	for k in range(0, len(Lines)):
		cut = Lines[k].split()
		current_line = k
		if(len(cut)>1):
			if(cut[0] == 'Begin' and cut[1] =='MD'):
				while(True):
					try:
						new_cut = Lines[current_line].split()
						if(len(new_cut) >1 and new_cut[0] =='Hirshfeld'):
							substrate_string = str(Lines[current_line+atom_number-5]).split()
							charge_string = str(Lines[current_line+atom_number+1]).split()
							if(len(atom_charge)<1):
								print('substrate: ' + substrate_string[2])
								print('plotting charge for: ' + charge_string[2])
								atom_letter = str(charge_string[2])

							current_charge = float(charge_string[1])

						elif(len(new_cut) >1 and new_cut[0]=='*' and new_cut[1]=='Maximum'):
							atom_charge.append(current_charge)
							break
						current_line +=1
						k +=1
					except:
						break
				#atom_charge.append(float(charge_string[1]))
			#current_line +=1

	atom_charge.remove(atom_charge[0])
	steps = list(range(0, len(atom_charge)))
	print(len(steps))
	fs_steps = [time_step*i for i in steps]
	return [fs_steps, atom_charge]


def plot_charge_single_outFile(path, atom_num):
	data = get_charge_single_outFile(path, atom_num)
	atom_charge = data[1]
	fig = px.line(x=data[0], y=data[1])
	fig.update_layout(
			title="Charge for atom: " + str(atom_num),
	    	xaxis_title="Time (fs)",
	    	yaxis_title="Charge (e)",
	    )

	fig.show()
	initial_avg = sum(atom_charge[0:20])/20
	print('intial charge: ' + str(initial_avg))
	max_charge = avg_max_charge(atom_charge)
	print("rolling avg max charge: " + str(max_charge))
	min_charge = avg_min_charge(atom_charge)
	print("rolling avg min charge: " + str(min_charge))
	avg = sum(atom_charge[-20:])/20
	delta_charge = avg - initial_avg
	print("final charge: " + str(avg))
	print("length of charge array: " +str(len(atom_charge)))
	print("delta_charge: " + str(delta_charge))
	#moving_average(atom_charge, "Charge (e)")

def plot_absolute_charge_single_outFile(path, atom_num):
	data = get_charge_single_outFile(path, atom_num)
	atom_charge = [abs(num) for num in data[1]]
	fig = px.line(x=data[0], y=atom_charge)
	fig.update_layout(
			title="Absolute Charge for atom: " + str(atom_num),
	    	xaxis_title="Time (fs)",
	    	yaxis_title="abs(Charge (e))",
	    )

	fig.show()
#---------------------------------------------------------------------------------------------------------------------------------#

#returns the sum of the charge change for 2 atoms and the fs time array
def get_charge_sum(path, atom1_number, atom2_number):
	file =open(path + '/out', 'r')
	Lines = file.readlines()
	current_line = 0
	counter_check =0
	atom1_charge = []
	atom2_charge = []
	time_step = 0.0242
	current_charge_atom1 = 0
	current_charge_atom2 = 0

	#only grabs the last hirshfeld atomic charge before next MD step, hence convoluted nested loops :)
	for k in range(0, len(Lines)):
		cut = Lines[k].split()
		current_line = k
		if(len(cut)>1):
			if(cut[0] == 'Begin' and cut[1] =='MD'):
				while(True):
					try:
						new_cut = Lines[current_line].split()
						if(len(new_cut) >1 and new_cut[0] =='Hirshfeld'):
							substrate_string = str(Lines[current_line+atom1_number-5]).split()
							charge1_string = str(Lines[current_line+atom1_number+1]).split()
							charge2_string = str(Lines[current_line+atom2_number+1]).split()
							
							if(len(atom1_charge)<1):
								print('substrate: ' + substrate_string[2])
								print('plotting charge for: ' + charge1_string[2] + ' and ' + charge2_string[2])
								atom1_letter = str(charge1_string[2])
								atom2_letter = str(charge2_string[2])
							
							current_charge_atom1 = float(charge1_string[1])
							current_charge_atom2 = float(charge2_string[1])

						elif(len(new_cut) >1 and new_cut[0]=='*' and new_cut[1]=='Maximum'): #only append last charge in scf loop; indicated by this string
							atom1_charge.append(current_charge_atom1)
							atom2_charge.append(current_charge_atom2)
							break
						current_line +=1
						k +=1
					except:
						print('error')
						break

	atom1_charge.remove(atom1_charge[0])
	atom2_charge.remove(atom2_charge[0])
	steps = list(range(0, len(atom1_charge)))
	print(len(steps))
	net_charge = np.add(atom1_charge, atom2_charge) #add two charge arrays together to get net change
	fs_steps = [time_step*i for i in steps]
	return [fs_steps, net_charge]

#returns charge array where every element is relative to initial charge 
def get_charge_sum_delta(path, atom1_number, atom2_number):
	data = get_charge_sum(path, atom1_number, atom2_number)
	charge_arr = data[1]
	initial_charge = charge_arr[0]
	charge_delta = [i - initial_charge for i in charge_arr]
	return [data[0], charge_delta]

#returns charge array where every element is relative to initial charge (used for single atom)
def get_charge_delta(path, atom_num):
	data = get_charge_single_outFile(path, atom_num)
	charge_arr = data[1]
	initial_charge = charge_arr[0]
	charge_delta = [i - initial_charge for i in charge_arr]
	return [data[0], charge_delta]

#Gets charge delta over range of atoms. for example from atom 56 to 60 (CH4 on 55 atom np)
def get_charge_delta_multiple_atoms(path, atom_start_number, atom_finish_number):
	num_atoms = atom_finish_number - atom_start_number +1
	atom_charges = []

	for atom in range(atom_start_number, atom_finish_number+1):
		data = get_charge_delta(path, atom)
		atom_charges.append(data[1])

	total_atoms_charge = []
	for i in range(len(atom_charges[0])):
		elem_sum = 0
		for j in range(num_atoms):
			elem_sum += atom_charges[j][i]
		
		total_atoms_charge.append(elem_sum)

	return [data[0], total_atoms_charge]
#---------------------------------------------------------------------------------------------------------------------#
def plot_charge_sum(path, atom1_number, atom2_number):
	data = get_charge_sum(path, atom1_number, atom2_number)
	fs_steps = data[0]
	net_charge = data[1]

	fig = px.line(x=fs_steps, y=net_charge)
	fig.update_layout(
			title='atom1 number: ' + str(atom1_number) + " atom2 number: " + str(atom2_number) + " Net Charge",
	    	xaxis_title="Time (fs)",
	    	yaxis_title="Charge (e)",
	    )
	#fig.add_vline(x=transition_state_time)
	fig.show()
	net_charge = net_charge.tolist()
	initial_avg = sum(net_charge[0:20])/20
	print('intial avg charge (over first 20 steps): ' + str(initial_avg))
	max_charge = avg_max_charge(net_charge)
	print("rolling avg max charge (avg between 10 points infront/after max value): " + str(max_charge))
	min_charge = avg_min_charge(net_charge)
	print("rolling avg min charge (avg between 10 points infront/after min value): " + str(min_charge))
	final_avg = sum(net_charge[-20:])/20
	print("final average charge: (over last 20 steps): " + str(final_avg))
	print(net_charge[-20:])
	print("length of charge array: " +str(len(net_charge)))
	delta_charge = final_avg - initial_avg
	biggest_oscillation = max(abs(min_charge - initial_avg), abs(max_charge - initial_avg))
	print("--------------------------")
	print('delta_charge: ' + str(delta_charge))
	print('biggest_oscillation: ' + str(biggest_oscillation))
#---------------------------------------------------------------------------------------------------------------------#

def get_bulk_charge_sum(path, atom1_start_number, atom2_finish_number):
	file =open(path + '/out', 'r')
	Lines = file.readlines()
	current_line = 0
	counter_check =0
	bulk_charge_at_time = []
	time_step = 0.0242
	md_step = 0
	end_md_step_count = 0

	#only grabs the last hirshfeld atomic charge before next MD step, hence convoluted nested loops :)
	for k in range(0, len(Lines)):
		cut = Lines[k].split()
		current_line = k
		if(len(cut)>1):
			if(cut[0] == 'Begin' and cut[1] =='MD'):
				md_step +=1
				while(True):
					try:
						new_cut = Lines[current_line].split()
						if(len(new_cut) >1 and new_cut[0] =='Hirshfeld' and new_cut[1]=='Net'):
							bulk_charge_sum = 0
							#print(new_cut)
							for i in range(atom1_start_number, atom2_finish_number +1):
								charge_string = str(Lines[current_line+i +1]).split()
								if len(bulk_charge_at_time)<1:
									print(charge_string)
								bulk_charge_sum += float(charge_string[1])
							current_line += atom2_finish_number

						elif(len(new_cut) >1 and new_cut[0]=='*' and new_cut[1]=='Maximum'):
							bulk_charge_at_time.append(bulk_charge_sum)
							break
						current_line +=1
						#k +=1
						k=current_line
					except Exception as e:
						print(e)
						break


	print(md_step)
	bulk_charge_at_time.remove(bulk_charge_at_time[0])
	steps = list(range(0, len(bulk_charge_at_time)))
	print(len(steps))
	fs_steps = [time_step*i for i in steps]
	return [fs_steps, bulk_charge_at_time]

#sum the charge change for 2 atoms and plot 
def plot_charge_sum_over_range(path, atom1_start_number, atom2_finish_number):
	file =open(path + '/out', 'r')
	Lines = file.readlines()
	current_line = 0
	counter_check =0
	bulk_charge_at_time = []
	time_step = 0.0242
	md_step = 0
	end_md_step_count = 0

	#only grabs the last hirshfeld atomic charge before next MD step, hence convoluted nested loops :)
	for k in range(0, len(Lines)):
		cut = Lines[k].split()
		current_line = k
		if(len(cut)>1):
			if(cut[0] == 'Begin' and cut[1] =='MD'):
				md_step +=1
				while(True):
					try:
						new_cut = Lines[current_line].split()
						if(len(new_cut) >1 and new_cut[0] =='Hirshfeld' and new_cut[1]=='Net'):
							bulk_charge_sum = 0
							#print(new_cut)
							for i in range(atom1_start_number, atom2_finish_number +1):
								charge_string = str(Lines[current_line+i +1]).split()
								#print(charge_string)
								bulk_charge_sum += float(charge_string[1])
							current_line += atom2_finish_number

						elif(len(new_cut) >1 and new_cut[0]=='*' and new_cut[1]=='Maximum'):
							bulk_charge_at_time.append(bulk_charge_sum)
							break
						current_line +=1
						#k +=1
						k=current_line
					except Exception as e:
						print(e)
						break


	print(md_step)
	steps = list(range(0, len(bulk_charge_at_time)))
	print(len(steps))
	fs_steps = [time_step*i for i in steps]
	fig = px.line(x=fs_steps, y=bulk_charge_at_time)
	fig.update_layout(
			title="total charge from atom " + str(atom1_start_number) + ' to ' + str(atom2_finish_number),
	    	xaxis_title="Time (fs)",
	    	yaxis_title="Charge (e)",
	    )
	fig.show()
	net_charge = bulk_charge_at_time
	initial_avg = sum(net_charge[0:20])/20
	print('intial avg charge (over first 20 steps): ' + str(initial_avg))
	max_charge = avg_max_charge(net_charge)
	print("rolling avg max charge (avg between 10 points infront/after max value): " + str(max_charge))
	min_charge = avg_min_charge(net_charge)
	print("rolling avg min charge (avg between 10 points infront/after min value): " + str(min_charge))
	final_avg = sum(net_charge[-20:])/20
	print("final average charge: (over last 20 steps): " + str(final_avg))
	print(net_charge[-20:])
	print("length of charge array: " +str(len(net_charge)))
	delta_charge = final_avg - initial_avg
	biggest_oscillation = max(abs(min_charge - initial_avg), abs(max_charge - initial_avg))
	print("--------------------------")
	print('delta_charge: ' + str(delta_charge))
	print('biggest_oscillation: ' + str(biggest_oscillation))
#-------------------------------------------------------------------------------------------------------------------#

#Get the charge sum at a given time step and avg over the last 20 points
def get_charge_sum_at_timestep(path, atom1_number, atom2_number, time_step_to_stop_at):
	file =open(path + '/out', 'r')
	Lines = file.readlines()
	current_line = 0
	counter_check =0
	atom1_charge = []
	atom2_charge = []
	time_step = 0.0242
	current_charge_atom1 = 0
	current_charge_atom2 = 0
	time_step_count = 0

	#only grabs the last hirshfeld atomic charge before next MD step, hence convoluted nested loops :)
	for k in range(0, len(Lines)):
		cut = Lines[k].split()
		current_line = k
		if(len(cut)>1):
			if(cut[0] == 'Begin' and cut[1] =='MD'):
				while(True):
					try:
						new_cut = Lines[current_line].split()
						if(len(new_cut) >1 and new_cut[0] =='Hirshfeld'):
							substrate_string = str(Lines[current_line+atom1_number-5]).split()
							charge1_string = str(Lines[current_line+atom1_number+1]).split()
							charge2_string = str(Lines[current_line+atom2_number+1]).split()
							#print(str(charge1_string) + " " + str(charge2_string))
							if(len(atom1_charge)<1):
								print('substrate: ' + substrate_string[2])
								print('plotting charge for: ' + charge1_string[2] + ' and ' + charge2_string[2])
								atom1_letter = str(charge1_string[2])
								atom2_letter = str(charge2_string[2])
							#if len(atom1_charge) < 100: print(charge1_string + " " + charge2_string)
							current_charge_atom1 = float(charge1_string[1])
							current_charge_atom2 = float(charge2_string[1])

						elif(len(new_cut) >1 and new_cut[0]=='*' and new_cut[1]=='Maximum'):
							atom1_charge.append(current_charge_atom1)
							atom2_charge.append(current_charge_atom2)
							time_step_count +=1
							break
						current_line +=1
						k +=1
					except:
						print('error')
						break
		if time_step_count==time_step_to_stop_at:
			print('stopping at: ' + str(time_step_count))
			break


	atom1_charge.remove(atom1_charge[0])
	atom2_charge.remove(atom2_charge[0])
	steps = list(range(0, len(atom1_charge)))
	print(len(steps))
	net_charge = np.add(atom1_charge, atom2_charge) #add two charge arrays together to get net change
	print('charge sum at this point: ' + str(net_charge[-1]))
	net_charge = net_charge.tolist()
	initial_avg = sum(net_charge[0:20])/20
	print('intial avg charge (over first 20 steps): ' + str(initial_avg))
	max_charge = avg_max_charge(net_charge)
	print("rolling avg max charge (avg between 10 points infront/after max value): " + str(max_charge))
	min_charge = avg_min_charge(net_charge)
	print("rolling avg min charge (avg between 10 points infront/after min value): " + str(min_charge))
	final_avg = sum(net_charge[-20:])/20
	print("final average charge: (over last 20 steps): " + str(final_avg))
	print(net_charge[-20:])
	print("length of charge array: " +str(len(net_charge)))
	delta_charge = final_avg - initial_avg
	biggest_oscillation = max(abs(min_charge - initial_avg), abs(max_charge - initial_avg))
	print("--------------------------")
	print('delta_charge: ' + str(delta_charge))
	print('biggest_oscillation: ' + str(biggest_oscillation))
	#return [fs_steps, net_charge]
#--------------------------------------------------------------------------------------------------------------------#


#returns the system spin (up - down) from atom_num1 to atom_num2
def get_spin_over_range(path, atom1_number, atom2_number):
	file =open(path + '/out', 'r')
	Lines = file.readlines()
	current_line = 0
	md_step =0
	up_spin_arr = []
	down_spin_arr = []
	system_spin_arr = []
	time_step = 0.0242

	for k in range(0, len(Lines)):
		cut = Lines[k].split()
		current_line = k
		if(len(cut)>1):
			if(cut[0] == 'Begin' and cut[1] =='MD'):
				while(True):
					try:
						new_cut = Lines[current_line].split()
						if(len(new_cut) >1 and new_cut[0] =='Spin-projected'):
							up_spin = 0
							down_spin =0
							for x in range(current_line + atom1_number +1, current_line + atom2_number + 2):
								spin_line = Lines[x].split()
								up_spin += float(spin_line[1])
								down_spin += float(spin_line[2])
								current_line +=1
							up_spin_arr.append(up_spin)
							down_spin_arr.append(down_spin)
							diff = up_spin - down_spin
							print('MD step = ' + str(md_step) + ' net spin = ' + str(diff))
							system_spin_arr.append(diff)
							break
					except:
						print('error')
						break
					current_line +=1
				md_step +=1
