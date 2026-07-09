#!/usr/bin/env python3
# coding: utf-8

#usage: ./make_cfinput-stable-v3.1.py [-h] [--verbose] --map MAP [--cmp CMP] [--alpha ALPHA] --npoints NPOINTS [--plot3D] [--plot2D] [--plotsphere] [--plotspherical] [--save] [--sampling SAMPLING]
#
#options:
#  -h, --help           show this help message and exit.
#  --verbose            Verbose mode. Default: NO.
#  --map MAP            3D map in npz format created with make_cumulative3D-stable-v3.1.py.
#  --cmp CMP            Color scheme. Default: YlOrBr.
#  --alpha ALPHA        Alpha value. Default: 0.1.
#  --npoints NPOINTS    Number of points in the weigth map.
#  --plot3D             Plots the 3D map. Default: NO.
#  --plot2D             Plots the 2D map. Default: NO.
#  --plotsphere         Plots the map on a sphere. Default: NO.
#  --plotspherical      Plots the 3D map in spherical coordinates. Default: NO.
#  --save               Saves the input file.
#  --sampling SAMPLING  Point sampling scheme (dist or max). Default: dist.

# Usage example: ./make_cfinput-stable-v3.1.py --verbose --map cumulative_reduced3Dmap-256nside-2000max-5.0step.npz --npoints 100000 --sampling dist --plotspherical --plotwedge 90

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from numba import njit, prange

import numpy as np
import healpy as hp

import argparse
import time
import sys
import os

import gc


def read_data():
	global number_pixel, distance_bin_number, distance_bins, distance_step, distance_max, dim
	global nside

	data = np.load(args.map)['cumulative_reduced3Dmap']
	
	#data = data[:48,:12]					#FOR TESTING 	
	#np.savetxt('test.csv',data,delimiter=',')		#FOR TESTING
	
	parameters = args.map.split('-')
	distance_max = int(parameters[2][:-3])
	distance_step = int(parameters[3][:-8])
	distance_bin_number = distance_max/distance_step
	
	if not distance_bin_number.is_integer():
		print('The number of bins does not seem to be an integer. Something really wrong in your map. I am quitting now!')
		sys.exit()
	else:
		distance_bin_number = int(distance_bin_number)
		
	nside = int(parameters[1][:-5])

	number_pixel = len(data)

	if distance_bin_number != len(data[1,:])-2:
		print('The number of bins I read from the map does not match the namefile. Something really wrong in your map. I am quitting now!')
		sys.exit()
	
	if nside != hp.npix2nside(number_pixel):
		print('The nside I read from the map does not match the namefile. Something really wrong in your map. I am quitting now!')
		sys.exit()	
	
	distance_bins = np.linspace(0,distance_max,distance_bin_number+1)

	dim = number_pixel * distance_bin_number
	#print(distance_bin_number);print(distance_bins)	#FOR TESTING

	if args.verbose:
		print('The 2D map has',str(number_pixel),'pixels and',str(distance_bin_number),'distance bins with size',str(distance_step),'Mpc from', \
		str(distance_bins[0]), 'Mpc to',str(distance_bins[-1]), 'Mpc for a total of',str(dim),'cells.')	

	return data


@njit(fastmath=True, error_model='numpy', boundscheck=True, parallel=True)
def select_data(data, npoints, sampling):

    flattened_data = data[:, 1:].T.flatten()

#######

    #distances = np.repeat(np.arange(distance_bin_number) * distance_step, number_pixel)

    n_pix = data.shape[0]
    n_dist = data.shape[1] - 1
    distances = np.repeat(np.arange(n_dist) * distance_step, n_pix)

#######

    total = np.sum(flattened_data)
    probabilities = flattened_data / total

    if sampling == 1:

        cdf = np.cumsum(probabilities)
        cdf /= cdf[-1]

        r = np.random.random(npoints)
        r = np.minimum(r, 1.0 - 1e-15)

        selected_points = np.searchsorted(cdf, r)

    else:

        npoints = min(npoints, flattened_data.shape[0])

        selected_points = np.argpartition(-flattened_data, npoints - 1)[:npoints]

    selected_weights = flattened_data[selected_points]
    selected_distances = distances[selected_points]

    return selected_points % n_pix, selected_distances, selected_weights

	
def dist2z(dist):
	H0 = 72  #km/s/MPc
	c = 3*10**5 #km/s
	return dist*H0/c

#@njit(fastmath=True,error_model='numpy',parallel=True) 
#def make_distance_bins(i,j,k):
#    dim = number_pixel * distance_bin_number
#    dist_i = np.empty(dim)
#    dist_j = np.empty(dim)
#    dist_k = np.empty(dim)
#    for n in prange(distance_bin_number):
#        dist_i[n*number_pixel:(n+1)*number_pixel] = i*n*distance_step
#        dist_j[n*number_pixel:(n+1)*number_pixel] = j*n*distance_step
#       dist_k[n*number_pixel:(n+1)*number_pixel] = k*n*distance_step
#   return dist_i,dist_j,dist_k	


def plot3Dmap(ras, decs, zs, weights):

	fig = plt.figure(figsize=(12, 12))
	ax3D = fig.add_subplot(projection='3d')		
	cm = plt.get_cmap(args.cmp)
		
	if args.verbose:
		print('Plotting the reduced 3D map (cumulative_3Dreducedmap-nside-max-step-points.png). It may take a while.')
		
	color_norm = mpl.colors.LogNorm(vmin=np.min(weights),vmax=np.max(weights))	
	colorbar = fig.colorbar(mpl.cm.ScalarMappable(norm=color_norm, cmap=args.cmp),ax=ax3D)
	colorbar.ax.set_ylim(np.min(weights),np.max(weights))
	
	alpha = 0.5*(weights-np.min(weights))/(np.max(weights)-np.min(weights))

	scatter = ax3D.scatter(ras, np.cos(decs), zs, c=weights, marker=',', cmap=cm, alpha=alpha)#, depthshade=False)
	#scatter = ax3D.tricontourf(ras, np.cos(decs), zs, marker='.', cmap=cm)

	rgba = [0,0,0,0]
#	plot_axis_max = args.max_distance 
	ax3D.set_xlim((0,2*np.pi))
	ax3D.set_ylim(-1,1)#(0,np.pi))
#	ax3D.set_zlim((-plot_axis_max,plot_axis_max))
	ax3D.set_xlabel('R.A. (rad)', fontsize=16)
	ax3D.set_ylabel('cos(dec)', fontsize=16)#('dec (rad)', fontsize=16)
		
	plt.savefig('cumulative_3Dreducedmap-'+str(nside)+'-'+str(distance_max)+'-'+str(distance_step)+'-'+str(args.npoints)+'.png')
	plt.close()

#@njit(fastmath=True,error_model='numpy',parallel=True) 
@njit(fastmath=True, error_model='numpy', boundscheck=True, parallel=False)
def spheric2cartesian(theta,phi):
	x = np.sin(theta)*np.cos(phi)
	y = np.sin(theta)*np.sin(phi)
	z = np.cos(theta)
	return x, y, z

def plot2Dmap(ras, decs, distances, weights):

	fig = plt.figure(figsize=(12, 12))
	ax2D = fig.add_subplot()		
	cm = plt.get_cmap(args.cmp)
		
	if args.verbose:
		print('Plotting the reduced 2D map (cumulative_2Dreducedmap-nside-max-step-npoints.png). It may take a while.')
		
	color_norm = mpl.colors.LogNorm(vmin=np.min(weights),vmax=np.max(weights))	
	colorbar = fig.colorbar(mpl.cm.ScalarMappable(norm=color_norm, cmap=args.cmp),ax=ax2D)
	colorbar.ax.set_ylim(np.min(weights),np.max(weights))
	
	alpha = 0.1*(distances-np.min(distances))/(np.max(distances)-np.min(distances))
	
	if np.isnan(alpha).any():
		alpha = 1
		print('Warning! Using alpha = 1 for the 2D map (all points are on a shell at fixed z).')
	
	size = 10**(2*(distances-np.max(distances))/(np.min(distances)-np.max(distances)))

	if np.isnan(size).any():
		size = 1
		print('Warning! Using size = 1 for the 2D map  (all points are on a shell at fixed z).')

	color = (weights-np.min(weights))/(np.max(weights)-np.min(weights))

	scatter = ax2D.scatter(ras, np.cos(decs), c=color, s=size, marker=',', cmap=cm, alpha=alpha,zorder=-100) #, depthshade=False)

	rgba = [0,0,0,0]
	ax2D.set_xlim((0,2*np.pi))
	ax2D.set_ylim(-1,1)#(0,np.pi))
	ax2D.set_xlabel('R.A. (rad)', fontsize=16)
	ax2D.set_ylabel('cos(dec)', fontsize=16)#('dec (rad)', fontsize=16)
		
	plt.savefig('cumulative_2Dreducedmap-'+str(nside)+'-'+str(distance_max)+'-'+str(distance_step)+'-'+str(args.npoints)+'.png')
	plt.close()

def plotspheremap(ras, decs, distances, weights):

	fig = plt.figure(figsize=(12, 12))
	axsphere = fig.add_subplot(projection='3d')		
	cm = plt.get_cmap(args.cmp)
		
	if args.verbose:
		print('Plotting the reduced map on a sphere (cumulative_spherereducedmap-nside-max-step-npoints.png). It may take a while.')
		
	color_norm = mpl.colors.LogNorm(vmin=np.min(weights),vmax=np.max(weights))	
	colorbar = fig.colorbar(mpl.cm.ScalarMappable(norm=color_norm, cmap=args.cmp),ax=axsphere)
	colorbar.ax.set_ylim(np.min(weights),np.max(weights))
	
	alpha = 0.2*((distances-np.max(distances))/(np.min(distances)-np.max(distances)))**2+0.2
	
	if np.isnan(alpha).any():
		alpha = 1
		print('Warning! Using alpha = 1 for the sphere map (all points are on a shell at fixed z).')
	
	size = 1 #10**(3*(zs-np.max(zs))/(np.min(zs)-np.max(zs)))

	if np.isnan(size).any():
		size = 1
		print('Warning! Using size = 1 for the sphere map  (all points are on a shell at fixed z).')

	color = (weights-np.min(weights))/(np.max(weights)-np.min(weights))

	x,y,z= spheric2cartesian(decs,ras)
	scatter = axsphere.scatter(x,y,z, c=color, s=size, marker=',', cmap=cm, alpha=alpha) #,zorder=-100*np.mean(zs))
	
	rgba = [0,0,0,0]
	axsphere.set_xlim(-1,1)
	axsphere.set_ylim(-1,1)
	#axsphere.set_xlabel('R.A. (rad)', fontsize=16)
	#axsphere.set_ylabel('cos(dec)', fontsize=16)#('dec (rad)', fontsize=16)
		
	plt.savefig('cumulative_spherereducedmap-'+str(nside)+'-'+str(distance_max)+'-'+str(distance_step)+'-'+str(args.npoints)+'.png')
	plt.close()

def plotsphericalmap(ras, decs, distances, weights):

	fig = plt.figure(figsize=(12, 12))
	axspherical = fig.add_subplot(projection='3d')		
	cm = plt.get_cmap(args.cmp)
		
	if args.verbose:
		print('Plotting the reduced spherical map (cumulative_sphericalreducedmap-nside-max-step-npoints.png). It may take a while.')
		
	color_norm = mpl.colors.LogNorm(vmin=np.min(weights),vmax=np.max(weights))	
	colorbar = fig.colorbar(mpl.cm.ScalarMappable(norm=color_norm, cmap=args.cmp),ax=axspherical)
	colorbar.ax.set_ylim(np.min(weights),np.max(weights))
	
	color = (weights-np.min(weights))/(np.max(weights)-np.min(weights))

	x,y,z = spheric2cartesian(decs,ras) * distances
	
	
	mask = ((x >= -distance_max/np.sqrt(2)) & (x <= distance_max/np.sqrt(2)) &(y >= -distance_max/np.sqrt(2)) & (y <= distance_max/np.sqrt(2)) \
		&(z >= -distance_max/np.sqrt(2)) & (z <= distance_max/np.sqrt(2)))
	x = x[mask]
	y = y[mask]
	z = z[mask]
	weights = weights[mask]
	color = color[mask]
	distances = distances[mask]

	alpha = 10**(2*(distances-np.min(distances))/(np.max(distances)-np.min(distances))-3)
		
	if np.isnan(alpha).any():
		alpha = 1
		print('Warning! Using alpha = 1 for the spherical map.')
		
	size = 10**(2*(distances-np.min(distances))/(np.max(distances)-np.min(distances)))+1
	scatter = axspherical.scatter(x,y,z, c=color, s=size, marker='.', cmap=cm, alpha=alpha) #,zorder=-100*np.mean(size))

#	size = 6
#	scatter = axspherical.scatter(x,y,z,c=color,s=size,marker='o',cmap=cm,alpha=0.25,linewidths=0)

	
	rgba = [0,0,0,0]
	axspherical.set_xlim(-distance_max/np.sqrt(2),distance_max/np.sqrt(2))
	axspherical.set_ylim(-distance_max/np.sqrt(2),distance_max/np.sqrt(2))
	axspherical.set_zlim(-distance_max/np.sqrt(2),distance_max/np.sqrt(2))
	
	#axsphere.set_xlabel('R.A. (rad)', fontsize=16)
	#axsphere.set_ylabel('cos(dec)', fontsize=16)#('dec (rad)', fontsize=16)
		
	plt.savefig('cumulative_sphericalreducedmap-'+str(nside)+'-'+str(distance_max)+'-'+str(distance_step)+'-'+str(args.npoints)+'.png')
	plt.close()
		
		
if __name__ == "__main__":

	global args

	parser = argparse.ArgumentParser(description='Program to make the CF input file and various plots of the 3D map.',prog='make_cfinput-stable-v3.1.py')
	parser.add_argument('--verbose', help='Verbose mode. Default: NO.', action='count')
	parser.add_argument('--map', help='3D map in npz format created with make_cumulative3D-stable-v3.1.py.', required=True)
	parser.add_argument('--cmp', help='Color scheme. Default: YlOrBr.', required=False, default='YlOrBr')
	parser.add_argument('--alpha', help='Alpha value. Default: 0.1.', required=False)
	parser.add_argument('--npoints', help='Number of points in the weight map.', required=True, type=int)
	parser.add_argument('--plot3D', help='Plots the 3D map. Default: NO.', action='count')
	parser.add_argument('--plot2D', help='Plots the 2D map. Default: NO.', action='count')
	parser.add_argument('--plotsphere', help='Plots the map on a sphere. Default: NO.', action='count')
	parser.add_argument('--plotspherical', help='Plots the 3D map in spherical coordinates. Default: NO.', action='count')
	parser.add_argument('--save', help='Saves the input file.', action='count')
	parser.add_argument('--sampling', help='Point sampling scheme (dist or max). Default: dist.', required=False, default='dist')
	args = parser.parse_args()
	
	check_filename = 'cf_input_' + args.map[:-4] + '-' + str(args.npoints) + '.txt'
			
	if os.path.exists(check_filename):
		print('Warning! The CF file (' + check_filename + ') already exists. Out of caution I will use it and not generate and not save a new CF file.')	    
		decs, ras, distances, weights = np.loadtxt(check_filename, unpack=True, skiprows=2,delimiter=',')
		check_filename_split = check_filename.split('-')
		nside = int(check_filename_split[1][:-5])
		distance_max = int(check_filename_split[2][:-3])
		distance_step = check_filename_split[3][:-4]
				
	else:
		if args.verbose:
			print('Reading the cumulative map.')
		
		data = read_data()
		
		if args.verbose:
			map_prob = np.round(100*np.sum(data[:,1:].T.flatten()),2)
			print(f'The map contains {map_prob}% of the total event probability.')	
			
		if args.npoints > dim:
			print('Warning! You are selecting a number of cells larger than their total number. If using --sampling max I will default to the latter. Please choose a value of --npoints smaller than', \
				str(args.npoints),'and rerun.')
			sys.exit()
		
		if args.verbose:
			print('Selecting ' + str(args.npoints) + ' pixels according to the weight probability distribution ('+str(np.round(args.npoints/dim*100,1))+'% of the total).')

		if args.sampling == 'dist':
			sampling = 1 
		elif args.sampling == 'max':
			sampling = 0
		else:
			print('I do not understand what sampling you want: dist or max? Try again. I am quitting now!')
			sys.exit()
		
		points, distances, weights = select_data(data,args.npoints,sampling)
	
		if args.verbose:
			print('Computing coordinates.')

		decs,ras = hp.pix2ang(nside, points, nest=True)
	
		if args.save:
			if args.verbose:
				print('Saving the input file.')
			if sampling == 1:
				header = 'Sample of ' + str(args.npoints) + ' probability pixels from ' + args.map + '\n' + 'dec,ra,z,weight'
			else:
				header = 'Highest ' + str(args.npoints) + ' probability pixels from ' + args.map + '\n' + 'dec,ra,z,weight'
		
			matrix_output = np.c_[decs,ras,distances,weights]
			matrix_output = matrix_output[matrix_output[:,3].argsort()][::-1]   #to sort in descending order
			np.savetxt('cf_input_'+args.map.split('.')[0]+'-'+str(args.npoints)+'.txt', matrix_output, fmt='%.6f,%.6f,%.6f,%.8e', delimiter=',', newline='\n', header=header, footer='', comments='# ', encoding=None)

		catch = np.round(100*np.sum(matrix_output[:,3]),2)

		if args.verbose:
			print(f'Your sampling is catching {catch}% of the total event probability.')

		if catch > 100:
			print('Warning! You are oversampling. The result may be unreliable. I suggest to decrease the value of --npoints.')	
		

		#print(ras);print(decs);print(points);print(distances);print(zs);print(weights)	#FOR TESTING

	if args.plot3D:
		plot3Dmap(ras, decs, distances, weights)
	
	if args.plot2D:
		plot2Dmap(ras, decs, distances, weights)
	
	if args.plotsphere:
		#zs = dist2z(distances)
		plotspheremap(ras, decs, distances, weights)

	if args.plotspherical:
		plotsphericalmap(ras, decs, distances, weights)


	sys.exit()
