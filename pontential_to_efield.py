import math
import numpy as np; np.random.seed(0)
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
from mpl_toolkits.mplot3d import axes3d
import pandas as pd
from statistics import mean
import os
from shutil import copyfile
import plotly.figure_factory as ff
from efield import *
import seaborn as sns; sns.set_theme()
from efield import *

###########-------VARIABLES TO CHANGE-----------###########
num_lines_in_x_direction = 101 #51 for normal calcs, 101 for higher res
num_slices_in_y_direction = 101
num_points_in_z_direction = 201 #101 for normal calcs, 201 for higher res
cell_length_bohr = 47.2431968887537 #Au55
cell_length_ang = 25 #Au55
#cell_length_bohr = 37.7945225 #Au19
#cell_length_ang = 20 #Au19
#cell_length_bohr = 66.1404144 #Au19 two NPS
#cell_length_ang = 35 #Au19 two NPS
num_vh_files = 60
heat_map_max_value = 0.3
heat_map_min_value = 0
heat_map_color = "Blues"
#heat_map_color= sns.color_palette("coolwarm")
abs_value_of_excess_field = True
calculate_excess = True
add_external_field = True #used if we need to add actual external field value 
num_time_steps = 3000
##########-------------------------------------############

def add_gs_na_potentials(file1, file2):
	f1 = open(file1, 'r')
	f2 = open(file2, 'r')
	f3 = open('combined', 'w')
	f1_lines = file1.readlines()
	f2_lines = file2.readlines()

	for i in range(len(f1_lines)):
		f1_cut = f1_lines[i].split()
		f2_cut = f2_lines[i].split()

		potential_sum = f1_cut[5] + f2_cut[5]
		coords_str = f1_cut[0] + ' ' + ' ' + f1_cut[1] + ' ' + f1_cut[2] + ' ' + f1_cut[3] + ' ' + f1_cut[4]
		f3.write(coords_str+ ' ' + str(potential_sum))

def plot_potential(filename, cell_length, direction_to_plot):
	data = get_potential(filename, cell_length, direction_to_plot)

	plt.plot(data[0], data[1])
	plt.legend()
	plt.xlabel("Distance in " + direction_to_plot + " [Angstrom]")
	plt.ylabel("Potential (Ry/e)")
	plt.show()

def get_potential(filename, cell_length, direction_to_plot):
	f = open(filename, 'r')
	Lines = f.readlines()
	y_array = [] #will hold difference from surface of NP to point in space where efield is calculated 
	z_array = [] #will hold distance along z direction of Au surface
	x_array = []
	potential_array = []

	current_y_value = 0.0
	temp_y_arr = []
	temp_potential_arr = []

	for i in range(len(Lines)):

		cut = Lines[i].split()
		x = cut[2]
		y = cut[3]
		z = cut[4]
		potential = float(cut[5])

		if float(y)==current_y_value:
			temp_y_arr.append(float(y))
			temp_potential_arr.append(potential)
		else:
			y_array.append(mean(temp_y_arr))
			potential_array.append(mean(temp_potential_arr))
			current_y_value = float(y)
			temp_y_arr.clear()
			temp_potential_arr.clear()
			count = 0

		x_array.append(float(x))
		#y_array.append(float(y))
		z_array.append(float(z))
		#potential_array.append(potential)

	y_array = [cell_length*i for i in y_array]
	z_array = [cell_length*i for i in z_array]
	x_array = [cell_length*i for i in x_array]

	if direction_to_plot == 'x':
		return [x_array, potential_array]
	elif direction_to_plot == 'y':
		return [y_array, potential_array]
	else:
		return [z_array, potential_array]

#cell_length in bohr radii, potential in Ry/e , array_to_return is direction (x,y,z) efield is calculated in
###Efield value is assigned at midpoint of 2 points given i.e. (x1+x2)/2, (y1+y2)/2, (z1+z2)/2 
def efield_from_potential(filename, cell_length, array_to_return):
	f = open(filename, 'r')
	Lines = f.readlines()
	y_array = [] #holds y coord of every efield point- calculated as midpoint between 2 points
	z_array = [] #holds z coord of every efield point- calculated as midpoint between 2 points
	x_array = [] #holds x coord of every efield point- calculated as midpoint between 2 points
	efield_array = []
	lattice_vector = 25 #angstrom
	i=0

	for i in range(len(Lines)-1):
		cut = Lines[i].split()
		x1 = cut[2]
		y1 = cut[3]
		z1 = cut[4]
		cut2 = Lines[i+1].split()
		x2 = cut2[2]
		y2 = cut2[3]
		z2 = cut2[4]
		fractional_coord1 = [float(x1), float(y1), float(z1)]
		fractional_coord2 = [float(x2), float(y2), float(z2)]
		coord1 = [i*cell_length for i in fractional_coord1] #convert fractional coordinates to Bohr 
		coord2 = [i*cell_length for i in fractional_coord2]
		potential_coord1 = float(cut[5])
		potential_coord2 = float(cut2[5])

		if float(z2) <float(z1):
			z_coord = 1
		else:
			x_coord = (float(x2) + float(x1))/2
			y_coord = (float(y2) + float(y1))/2
			z_coord = (float(z2) + float(z1))/2
		y_array.append(y_coord)
		z_array.append(z_coord)
		x_array.append(x_coord)

		#a,b,c are 3 values which are the differnece in position x2-x1, y2-y1, z2-z1, respectively
		a = float(coord2[0]) - float(coord1[0])
		b = float(coord2[1]) - float(coord1[1])
		c = float(coord2[2]) - float(coord1[2])
		length = math.sqrt((a*a)  +(b*b) + (c*c)) #now calculate distance between our 2 points (x1,y1,z1) and (x2,y2,z2)

		efield = (potential_coord2 - potential_coord1)/length #units of Ry/Bohr/e
		efield_array.append(efield)
		#print('efield between [' +x1 +',' + y1 +',' +z1 +'] and ['+x2 +',' + y2 +',' +z2 +'] = ' + str(efield))
		#print('efield defined at coordinates: ' + str(x_coord) + ', ' + str(y_coord) + ', ' + str(z_coord))
		i +=1
	
	y_array = [lattice_vector*i for i in y_array] #scale arrays to units of angstrom 
	z_array = [lattice_vector*i for i in z_array]
	x_array = [lattice_vector*i for i in x_array]

	if array_to_return == 'x':
		return [x_array, np.array(efield_array)]
	elif array_to_return =='y':
		return [y_array, np.array(efield_array)]
	else:
		return [z_array, np.array(efield_array)]

def overlay_plots(is_matplotlib):
	filenames = ['y_out_30_0-66x_0-49z', 'y_out_35']
	time_points = ['x=0.66, z=0.49', 'x=0.49, z=0.49']
	distances = []
	efields = []
	direction = 'y'

	for file in filenames:
		data = efield_from_potential(file, cell_length_bohr, direction)
		distances.append(data[0])
		efields.append(data[1])

	if is_matplotlib:
		for i in range(len(filenames)):
			plt.plot(distances[i], efields[i], label = time_points[i])
		plt.legend()
		plt.xlabel("distance along unit cell (ang)")
		plt.ylabel("calculated efield strength (Ry/Bohr/e)")
		plt.show()
	else:
		fig = go.Figure()
		for i in range(len(filenames)):
			fig.add_trace(go.Scatter(x=distances[i], y=efields[i],
							mode='lines', name=time_points[i]))

		fig.update_layout(title='Calculated E-field Above Au', xaxis_title='Distance above Au Surface [Angstrom]', 
			yaxis_title='Calculated E-field Strength (Ry/Bohr/e)')
		fig.update_xaxes(range=[1, 6])
		fig.show()

#plots the difference between two efields (i.e. 1 in the ground state and 1 with external field present) **subtracts efield of 2nd file from the first 
def difference(file1, file2):
	#file1 = 'Au_noEfield/y_out_00'
	#file2 = 'Au55_0.1_fixed/y_out_45'
	efield1 = []
	efield2 = []
	distances = []
	direction = 'y'

	data1 = efield_from_potential(file1, cell_length_bohr, direction)
	distances = data1[0]
	efield1 = data1[1]
	data2 = efield_from_potential(file2, cell_length_bohr, direction)
	efield2 = data2[1]

	diff = np.subtract(efield1, efield2)
	#return [distances, diff]
	plt.plot(distances, diff)
	#plt.axvline(x=6.81, c='red') #y 
	#plt.axvline(x=15.38, c='red') #y
	#plt.axvline(x=8.13, c='red') #z 
	#plt.axvline(x=16.45, c='red') #z
	plt.legend()
	plt.xlabel("Distance along unit cell (Angstrom)")
	plt.ylabel("Difference in efield strength (Ry/Bohr/e)")
	plt.yticks(np.arange(-0.7, 0.7, 0.1))
	plt.show()

def average_over_line():
	o2_heights = [2.5, 3, 3.4, 3.8, 4] 
	radius_below = 2.4 #region to average over above and below O2 height in angstrom
	radius_above = 0
	file = 'y_out_50'
	direction = 'y'
	y_surface = 15.38 #angstrom

	data = efield_from_potential(file, cell_length_bohr, direction)
	distance_array = data[0] #distances in angstrom
	efield = data[1]
	print(len(efield))

	for i in range(len(o2_heights)):
		lower_bound = y_surface + o2_heights[i] - radius_below #units of angstrom
		upper_bound = y_surface + o2_heights[i] + radius_above #angstrom
		closest_lower_bound = get_closest(lower_bound, distance_array) #index of the closest lower bound
		closest_upper_bound = get_closest(upper_bound, distance_array) #index of the closest upper bound
		print("averaging between: " + str(distance_array[closest_lower_bound]) + " and "+ str(distance_array[closest_upper_bound]))
		average_efield = mean(efield[closest_lower_bound:closest_upper_bound])
		print(average_efield)

#returns the index of given array which contains value closest to given point
def get_closest(point, array):
	distance_from_array = [abs(i-point) for i in array]
	min_value = min(distance_from_array)
	return distance_from_array.index(min_value)

#calculate efield along a line (direction determined by points used in grid2val), repeat for all lines on plane
#returned array is difference between case with external field and no external field along a line, for every line on plane
def average_over_plane(external_field_file_name, zero_field_file_name, calculate_excess_bool, direction):
	external_field_file = open(external_field_file_name, "r")
	external_field_file_lines = external_field_file.readlines()
	zero_field_file = open(zero_field_file_name, "r")
	zero_field_file_lines = zero_field_file.readlines()

	direction = 'z'
	Efield_array = []
	current_line_position = 0

	print('avg over plane, before while')

	while current_line_position <len(external_field_file_lines):
		temporary_external_field_file = open("temporary_ext_field", "w")
		temporary_zero_field_file = open("temporary_zero_field", "w")

		for k in range(num_points_in_z_direction):
			external_field_line = external_field_file_lines[current_line_position]
			zero_field_line = zero_field_file_lines[current_line_position]
			temporary_external_field_file.write(external_field_line)
			temporary_zero_field_file.write(zero_field_line)
			current_line_position +=1

		temporary_external_field_file.close()
		temporary_zero_field_file.close()

		external_field_data = efield_from_potential("temporary_ext_field", cell_length_bohr, direction)
		zero_external_field_data = efield_from_potential("temporary_zero_field", cell_length_bohr, direction)

		if calculate_excess_bool:
			difference_in_fields = []
			for i in range(len(zero_external_field_data[1])):
				difference_in_fields.append(external_field_data[1][i] - zero_external_field_data[1][i])

			abs_difference = [abs(ele) for ele in difference_in_fields]

			if abs_value_of_excess_field:
				Efield_array.append(abs_difference)
			else:
				Efield_array.append(difference_in_fields) #if we want to see negative efield values too
		else:
			abs_efield = [abs(ele) for ele in external_field_data[1]]
			if abs_value_of_excess_field:
				Efield_array.append(abs_efield)
			else:
				Efield_array.append(external_field_data[1])

	print('finished avg over plane')
	return Efield_array	

def heat_map_over_time(x_axis_name, y_axis_name):
	direction = 'z'
	start_point = 0
	end_point = 60
	#end_point = num_vh_files 
	max_efield_at_time = []
	avg_efield_at_time = []
	for i in range(start_point, end_point+1):
		data = average_over_plane('allPoints_out_' + str(i), 'allPoints_out_noField_0', calculate_excess, direction)
		print(len(data[0]))
		print(len(data))
		fs_time = round((i/num_vh_files)*num_time_steps*.0242, 1)

		for z_array in data: #have to insert 0 into first position of every array since 1st efield point is
			z_array.insert(0, 0)

		if add_external_field:
			efield_value = calculate_efield_at(fs_time)
			print('value of external field at time: ' + str(fs_time) + ' is: ' + str(efield_value))
			for x in range(len(data)):
				for y in range(len(data[x])):
					data[x][y] += abs(efield_value)

		max_values = []
		avg_values = []
		for arr in data:
			max_values.append(max(arr))
			avg_values.append(mean(arr))

		max_efield_at_time.append(max(max_values))
		avg_efield_at_time.append(mean(avg_values))
		print('t= ' + str(fs_time) + ' max efield: ' + str(max(max_values)))
		print('t= ' + str(fs_time) + ' avg efield: ' + str(mean(avg_values)))

		##convert to data frame##
		num_columns = len(data[0])
		num_rows = len(data)
		print('columns: ' + str(num_columns))
		print('rows: ' + str(num_rows))
		list_of_columns = [round((j/(num_columns))*cell_length_ang, 0) for j in range(num_columns)]
		print('last element in columns: ' + str(list_of_columns[-1]))
		list_of_rows = [round((j/(num_rows))*cell_length_ang, 0) for j in range(num_rows)]
		#list_of_rows = list_of_rows.reverse()
		print('last element in rows: ' + str(list_of_rows[-1]))
		df = pd.DataFrame(data, columns = [list_of_columns], index=[list_of_rows])

		#-------------------make heat map from dataframe--------------------------------------------
		ax = sns.heatmap(df, xticklabels = 15, yticklabels =10, vmin=heat_map_min_value, vmax=heat_map_max_value, cmap=heat_map_color)
		ax.set_xlabel(x_axis_name + " ($\AA$)")
		ax.set_ylabel(y_axis_name + " ($\AA$)")
		if calculate_excess:
			#ax.set_title("Excess E-Field at t=" + str(fs_time) + " and x= 12.5 ($\AA$)") #for fixed x value
			ax.set_title("Excess E-Field at t=" + str(fs_time) + " and z= 15.2 ($\AA$)") #for fixed height
		else:
			ax.set_title("Total E-Field at t=" + str(fs_time))
		plt.savefig('VH_' + str(i) + '.png', dpi=299)
		plt.clf()
	
	print('max efields:')
	print(max_efield_at_time)
	print('avg efields: ')
	print(avg_efield_at_time)
		#------------------------------------------------------------------------------------------
def heat_map_brain_slice(external_field_file, zero_field_file, time, height, direction):
	y_height = round((height/(num_slices_in_y_direction-1))*cell_length_ang, 1)
	data = average_over_plane(external_field_file, zero_field_file, calculate_excess, direction)
	print('-----Calculated Efield-------------')
	if height<5:
		print(data)
	print('--------------------------------')
	print(len(data[0]))
	print(len(data))
	fs_time = round((time/num_vh_files)*num_time_steps*.0242, 1)
	print('fs time: ' + str(fs_time))

	for z_array in data: #have to insert 0 into first position of every array since 1st efield point is
		z_array.insert(0, 0)

	if add_external_field:
			efield_value = calculate_efield_at(fs_time)
			print('value of external field at time: ' + str(fs_time) + ' is: ' + str(efield_value))
			for x in range(len(data)):
				for y in range(len(data[x])):
					data[x][y] += abs(efield_value)
					
	##convert to data frame##
	num_columns = len(data[0])
	num_rows = len(data)
	print('columns: ' + str(num_columns))
	print('rows: ' + str(num_rows))

	list_of_columns = [round((i/(num_columns))*cell_length_ang, 0) for i in range(num_columns)]
	print('last element in columns: ' + str(list_of_columns[-1]))
	list_of_rows = [round((i/(num_rows))*cell_length_ang, 0) for i in range(num_rows)]
	#list_of_rows = list_of_rows.reverse()
	print('last element in rows: ' + str(list_of_rows[-1]))
	df = pd.DataFrame(data, columns = [list_of_columns], index=[list_of_rows])

	#-------------------make heat map from dataframe--------------------------------------------
	ax = sns.heatmap(df, xticklabels = 15, yticklabels =10, vmin=heat_map_min_value, vmax=heat_map_max_value, cmap=heat_map_color)
	ax.set_xlabel("Z-Coordinate ($\AA$)")
	ax.set_ylabel("X-Coordinate ($\AA$)")
	if calculate_excess:
		ax.set_title("Excess E-Field at t=" + str(fs_time) + " and y= " + str(y_height) + " ($\AA$)")
	else:
		ax.set_title("Total E-Field at t=" + str(fs_time) + " and y= " + str(y_height) + " ($\AA$)")
	plt.savefig(str(height) + '.png', dpi=299)
	plt.clf()
	#------------------------------------------------------------------------------------------
	return [np.amax(data), y_height] #returns max value of calculated efield at given plane height 

def brain_slice(wfz_time_point):
	print('calculate_excess: ' + str(calculate_excess))
	print('calculate_abs value: ' + str(abs_value_of_excess_field))
	direction = 'z'
	num_points_per_plane = num_lines_in_x_direction*num_points_in_z_direction
	max_efield_at_height_list = []
	mean_efield_at_height_list = []
	list_of_heights = []

	external_field_file = open("allPoints_out_" + str(wfz_time_point), 'r')
	zero_field_file = open("allPoints_out_noField_" + str(wfz_time_point), 'r')
	external_field_file_lines = external_field_file.readlines()
	zero_field_file_lines = zero_field_file.readlines()
	y_height = 0 #

	print(len(external_field_file_lines))
	current_line_position = 0

	while current_line_position <len(external_field_file_lines):
		temporary_external_field_file = open("temporary_ext_field_brainSlice", "w")
		temporary_zero_field_file = open("temporary_zero_field_brainSlice", "w")
		for i in range(num_points_per_plane):
			external_field_line = external_field_file_lines[current_line_position]
			zero_field_line = zero_field_file_lines[current_line_position]
			temporary_external_field_file.write(external_field_line)
			temporary_zero_field_file.write(zero_field_line)
			current_line_position +=1

		temporary_external_field_file.close()
		temporary_zero_field_file.close()
		data = heat_map_brain_slice("temporary_ext_field_brainSlice", "temporary_zero_field_brainSlice", wfz_time_point, y_height, direction)
		#max_efield_at_height_list.append(max(data[0]))
		max_efield_at_height_list.append(data[0])
		print('max value at height: ' +str(data[1]) + ' Ang is: ' + str(data[0]))
		#mean_efield_at_height_list.append(mean(data[0]))
		list_of_heights.append(data[1])
		#break #DELETE THIS LINE AFTER TESTING
		y_height +=1

	#plt.plot(list_of_heights, max_efield_at_height_list)
	#plt.xlabel("Height in Y direction (Angstrom")
	#plt.ylabel("Max Efield Ry/Bohr/e")
	#plt.ylim([0, 0.5])
	#plt.savefig(str(wfz_time_point) + '.png', dpi=299)
	#plt.clf()

	#ads_vs_efield(max_efield_at_height_list, list_of_heights)
	

def max_field_plots():

	for i in range(20, 61):
		print("-----------------RUNING NUMBER: " + str(i) + " -------------------")
		brain_slice(i)

#max_field_plots()

#heat_map_over_time('Z-Coordinate', 'Y-Coordinate') #for fixed x
heat_map_over_time('Y-Coordinate', 'X-Coordinate') #for fixed height 
