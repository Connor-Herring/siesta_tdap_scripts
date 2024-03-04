import numpy as np
import numpy.fft as npf 
import matplotlib.pyplot as plt
import math
from statistics import mean
import plotly.express as px
import plotly.graph_objects as go

def plot_total_forces(path):
	f = open(path + '/out', 'r')
	avg_x_forces = []
	avg_y_forces = []
	avg_z_forces = []
	time = []
	have_num_atoms = False
	num_atoms =0
	steps = 0
	current_line =1

	Lines = f.readlines()
	for line in Lines:
		cut = line.split()

		if(len(cut)>1):

			if(have_num_atoms ==False and cut[0] =='NumberOfAtoms'):
				num_atoms = int(cut[1])
				have_num_atoms = True

			if(cut[0] == 'siesta:' and cut[1] == 'Atomic' and cut[2] =='forces'):
				x_temp = []
				y_temp = []
				z_temp = []
				for i in range(0, num_atoms):
					forces = Lines[current_line + i]
					forces_cut = forces.split()
					x_temp.append(abs(float(forces_cut[1])))
					y_temp.append(abs(float(forces_cut[2])))
					z_temp.append(abs(float(forces_cut[3])))

				avg_x_forces.append(mean(x_temp))
				avg_y_forces.append(mean(y_temp))
				avg_z_forces.append(mean(z_temp))
				time.append(float(steps*.0242))
				steps +=1

		current_line +=1

	plt.xlabel('t (fs)')
	plt.ylabel('F_tot (ev/Ang)')
	plt.plot(time, avg_x_forces)
	plt.plot(time, avg_y_forces)
	plt.plot(time, avg_z_forces)
	plt.show()

def plot_single_force(path, num_files, atom_number):
	count = 0
	fill_num = 0 #used for correct number of 0s in out filenames 
	atom_force_x = []
	atom_force_y = []
	atom_force_z = []

	for x in range(0, num_files+1):
		if(count<10):
			filename = path '/out.000' + str(fill_num) + str(count)
			f = open(filename, 'r')
			force_data = get_force(f, atom_number)
			atom_force_x.append(float(force_data[2]))
			atom_force_y.append(float(force_data[3]))
			atom_force_z.append(float(force_data[4]))
			#print('plotting force for: ' + force_data[0])
		else: 
			filename = path + '/out.000' + str(count)
			f = open(filename, 'r')
			force_data = get_force(f, atom_number)
			atom_force_x.append(float(force_data[2]))
			atom_force_y.append(float(force_data[3]))
			atom_force_z.append(float(force_data[4]))
			#print('plotting force for: ' + force_data[0])
		count +=1

	steps = list(range(0, len(atom_force_x)))*1.2 #72 fs/60 WFZ
	plt.plot(steps, atom_force_x, label='x')
	plt.plot(steps, atom_force_y, label='y')
	plt.plot(steps, atom_force_z, label='z')
	plt.legend()
	plt.xlabel('Time (fs)')
	plt.ylabel('force (ev/Ang)')
	plt.show()

def plot_force(path, atom_number):
	f = path
	data = get_force(f, atom_number)

	x_fig = px.line(x=data[3], y=data[0])
	x_fig.update_layout(
			title=" X-Force vs time",
	    	xaxis_title="Time (fs)",
	    	yaxis_title="Force (ev/Ang)",
	    )
	x_fig.show()

	max_x_force = max(data[0][0:-1])
	print('Max x force: ' + str(max_x_force))
	min_x_force = min(data[0][0:-1])
	print('Min x force: ' + str(min_x_force))

	y_fig = px.line(x=data[3], y=data[1])
	y_fig.update_layout(
			title=" Y-Force vs time",
	    	xaxis_title="Time (fs)",
	    	yaxis_title="Force (ev/Ang)",
	    )
	y_fig.show()

	max_y_force = max(data[1][0:-1])
	print('Max y force: ' + str(max_y_force))
	min_y_force = min(data[1][0:-1])
	print('Min y force: ' + str(min_y_force))

	z_fig = px.line(x=data[3], y=data[2])
	z_fig.update_layout(
			title=" Z-Force vs time",
	    	xaxis_title="Time (fs)",
	    	yaxis_title="Force (ev/Ang)",
	    )
	z_fig.show()

	max_z_force = max(data[2][0:-1])
	print('Max z force: ' + str(max_z_force))
	min_z_force = min(data[2][0:-1])
	print('Min z force: ' + str(min_z_force))

#searches through out file to find the force on a given atom, returns that line as a string array
def get_force(path, atom_number):
	file = open(path + '/out')
	Lines = file.readlines()
	current_line = 0
	atom_force_x =[]
	atom_force_y =[]
	atom_force_z =[]
	time_step = .0242
	test_count = 0

	for line in Lines:
		cut = line.split()
		if(len(cut)>1):
			if(cut[0] == 'siesta:' and cut[1] == 'Atomic' and cut[2] =='forces'):
				force_string = str(Lines[current_line+atom_number]).split()
				#if test_count <10:
					#print(force_string)
					#test_count +=1
				a = float(force_string[1])
				b = float(force_string[2])
				c = float(force_string[3])
				atom_force_x.append(a)
				atom_force_y.append(b)
				atom_force_z.append(c)
				
		current_line +=1

	steps = list(range(0, len(atom_force_x)))
	fs_time = [time_step*i for i in steps]

	return [fs_time, atom_force_x, atom_force_y, atom_force_z]

#return just the average force over a time range, time given in fs
def get_avg_force(path, atom_number, start_time, end_time):
	data = get_force(path, atom_number)

	start_index = int(start_time/.0242)
	end_index = (int(end_time/.0242))

	avg_x_force = mean(data[0][start_index:end_index])
	avg_y_force = mean(data[1][start_index:end_index])
	avg_z_force = mean(data[2][start_index:end_index])

	print('Avg x force over time range: ' + str(avg_x_force))
	print('Avg y force over time range: ' + str(avg_y_force))
	print('Avg z force over time range: ' + str(avg_z_force))

#get force over range of atoms: sum their repective x,y,z components
def get_force_sum(path, atom_start_number, atom_end_number):
	file = open(path + '/out')
	Lines = file.readlines()
	current_line = 0
	atom_force_x =[]
	atom_force_y =[]
	atom_force_z =[]
	time_step = .0242

	for line in Lines:
		cut = line.split()
		if(len(cut)>1):
			if(cut[0] == 'siesta:' and cut[1] == 'Atomic' and cut[2] =='forces'):
				temp_x_force = []
				temp_y_force = []
				temp_z_force = []
				for i in range(current_line + atom_start_number, current_line + atom_end_number +1):
					force_string = str(Lines[i]).split()
					if len(atom_force_x) <1:
						print(force_string)
					a = float(force_string[1])
					b = float(force_string[2])
					c = float(force_string[3])
					temp_x_force.append(a)
					temp_y_force.append(b)
					temp_z_force.append(c)

				atom_force_x.append(sum(temp_x_force))
				atom_force_y.append(sum(temp_y_force))
				atom_force_z.append(sum(temp_z_force))
				
		current_line +=1

	steps = list(range(0, len(atom_force_x)))
	fs_time = [time_step*i for i in steps]

	return [atom_force_x, atom_force_y, atom_force_z, fs_time]

def plot_force_over_range(path, atom_start_number, atom_end_number):
	f = path + '/out'
	data = get_force_sum(f, atom_start_number, atom_end_number)

	x_fig = px.line(x=data[3], y=data[0])
	x_fig.update_layout(
			title=" X-Force vs time",
	    	xaxis_title="Time (fs)",
	    	yaxis_title="Force (ev/Ang)",
	    )
	x_fig.show()

	y_fig = px.line(x=data[3], y=data[1])
	y_fig.update_layout(
			title=" Y-Force vs time",
	    	xaxis_title="Time (fs)",
	    	yaxis_title="Force (ev/Ang)",
	    )
	y_fig.show()

	z_fig = px.line(x=data[3], y=data[2])
	z_fig.update_layout(
			title=" Z-Force vs time",
	    	xaxis_title="Time (fs)",
	    	yaxis_title="Force (ev/Ang)",
	    )
	z_fig.show()
