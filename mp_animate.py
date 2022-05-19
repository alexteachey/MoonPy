import matplotlib.pyplot as plt 
import numpy as np
from astropy.io import fits 
from matplotlib import animation
from matplotlib import cm 

"""
Matplotlib Animation Example

author: Jake Vanderplas
email: vanderplas@astro.washington.edu
website: http://jakevdp.github.com
license: BSD
Please feel free to use and modify this, but keep the above information. Thanks!
"""


# First set up the figure, the axis, and the plot element we want to animate
"""
fig = plt.figure()
ax = plt.axes(xlim=(0, 2), ylim=(-2, 2))
line, = ax.plot([], [], lw=2)


# initialization function: plot the background of each frame
def init():
	line.set_data([], [])
	return line,

# animation function.  This is called sequentially
def animate(i):
	x = np.linspace(0, 2, 1000)
	y = np.sin(2 * np.pi * (x - 0.01 * i))
	line.set_data(x, y)
	return line,

"""

def pix2coord(pixel, ref_pixel, ref_coord, pixel_size_deg):
	#print('calling pix2coord.')
	#### takes a given pixel value and gives you the coordinate value.
	delta_pixel = pixel - ref_pixel  ### degrees
	delta_coord = delta_pixel * pixel_size_deg 
	output_coord = delta_coord + ref_coord
	return output_coord 

def coord2pix(coord, ref_pixel, ref_coord, pixel_size_deg):
	#print('calling coord2pix.')
	#### takes in a coordinate and gives you the pixel value
	delta_coord = coord - ref_coord 
	delta_pixel = delta_coord / pixel_size_deg 
	output_pixel = delta_pixel + ref_pixel 
	return output_pixel 




# initialization function: plot the background of each frame
"""
def init():
	#im.set_data(np.random.random((5,5))) ### ORIGINAL
	im.set_data(flux[0]) 
	fig.colorbar(im)
	ax.invert_xaxis() #### because RA increases to the left (EAST)	
	return [im]


# animation function.  This is called sequentially
def animate(i, im, fig, ax, times, flux):

	if i == 0:
		im.set_data(flux[0]) 
		fig.colorbar(im)
		ax.invert_xaxis() #### because RA increases to the left (EAST)	

	#a=im.get_array() #### ORIGINAL
	#a=a*np.exp(-0.001*i)    # exponential decay of the values ### ORIGINAL 
	next_frame = flux[i]
	#im.set_array(a) #### ORIGINAL
	ax.set_title('BKJD = '+str(round(times[i],2)))
	#ax.scatter(refx_pix, refy_pix, color='red', marker='*', s=100)
	ax.scatter(objx_pix, objy_pix, color='red', marker='*', s=100)
	try:
		ax.set_xticks(ticks=xticks, labels=xtick_coords)
		ax.set_yticks(ticks=yticks, labels=ytick_coords)
	except:
		ax.set_xticks(ticks=xticks)
		ax.set_yticks(ticks=yticks)
		ax.set_xticklabels(labels=xtick_coords)
		ax.set_yticklabels(labels=ytick_coords)		


	ax.set_xlabel('RA [deg]')
	ax.set_ylabel('Dec [deg]')
	im.set_array(next_frame)
	return [im]


"""



def animate_TESS_FFI(filedir=None, filename=None, filepath=None, save_animation=False, return_arrays=True, normalize_flux=True):
	#### test file 
	#filedir = '/Users/hal9000/Documents/Projects/Unistellar_transit_obs/TIC121114997_tesscut'
	#testfile = 'tess-s0014-2-3_287.154269_41.565792_10x15_astrocut.fits'
	#testfile = 'tess-s0040-2-4_287.154269_41.565792_10x15_astrocut.fits'
	#testfits = fits.open(filedir+'/'+testfile)

	if type(filepath) == type(None):
		filepath = filedir+'/'+filename 

	else:
		### filepath is specified, so you want to separate out the filedir and filename
		slash_idxs = []
		for str_idx, str_char in enumerate(filepath):
			if str_char == '/':
				slash_idxs.append(str_idx)
		filedir = filepath[:slash_idxs[-1]] #### going up the filepath and stopping before the last slash 
		filename = filepath[slash_idxs[-1]+1:] #### starting after the slash and going to the end

	fitsfile = fits.open(filepath)
	datapage = fitsfile[1].data
	pixelpage = fitsfile[2].header
	flux = datapage['FLUX']
	times = datapage['TIME']


	tess_pixel_size_arcseconds = 21
	tess_pixel_size_degrees = tess_pixel_size_arcseconds / (60 * 60) #### 60 arcseconds per arcminute, 60 arcminutes per degree

	refx_pix, refy_pix = pixelpage['CRPIX1'], pixelpage['CRPIX2'] #### pixel value of the reference point (should be the target?)
	refx_RA, refy_Dec = pixelpage['CRVAL1'], pixelpage['CRVAL2'] 
	obj_RA, obj_Dec = pixelpage['RA_OBJ'], pixelpage['DEC_OBJ']
	objx_pix = ((obj_RA - refx_RA) / tess_pixel_size_degrees) + refx_pix 
	objy_pix = ((obj_Dec - refy_Dec) / tess_pixel_size_degrees) + refy_pix

	objx_pix = coord2pix(coord=obj_RA, ref_pixel=refx_pix, ref_coord=refx_RA, pixel_size_deg=tess_pixel_size_degrees)
	objy_pix = coord2pix(coord=obj_Dec, ref_pixel=refy_pix, ref_coord=refy_Dec, pixel_size_deg=tess_pixel_size_degrees)

	print("ref RA, Dec (x,y) = "+str(refx_RA)+', '+str(refy_Dec)+' ('+str(refx_pix)+', '+str(refy_pix)+')')
	print('obj RA, Dec (x,y) = '+str(obj_RA)+', '+str(obj_Dec)+' ('+str(objx_pix)+', '+str(objy_pix)+')')


	#### remember that numpy arrays are indexed as array[row][column]
	### row is the y-coordinate! column is the x-coordinate!

	nrows = flux.shape[1] #### number of boxes along the y-axis!
	ncols = flux.shape[2] #### number of boxes along x-axis!

	xlim = (-0.5, ncols-0.5) #### (0,0) is in the CENTER of the 0,0th pixel, so you need to set the limits wider
	ylim = (-0.5, nrows-0.5) 

	xticks = np.linspace(xlim[0], xlim[1], 5)
	yticks = np.linspace(ylim[0], ylim[1], 5)

	xtick_coords = pix2coord(pixel=xticks, ref_pixel=refx_pix, ref_coord=refx_RA, pixel_size_deg=tess_pixel_size_degrees)
	ytick_coords = pix2coord(pixel=yticks, ref_pixel=refy_pix, ref_coord=refy_Dec, pixel_size_deg=tess_pixel_size_degrees)

	xtick_coords = np.around(a=xtick_coords, decimals=3)
	ytick_coords = np.around(a=ytick_coords, decimals=3)

	fig = plt.figure(figsize=(8,8))
	ax = plt.axes(xlim=xlim, ylim=ylim)

	#### alex's modification
	a = flux[0] #### start with the first page
	#a = flux[0] / np.nanmedian(flux, axis=0)
	im = plt.imshow(a, vmin=np.nanmin(flux), vmax=np.nanmax(flux), interpolation='none', origin='lower')



	# initialization function: plot the background of each frame
	def init():
		#im.set_data(np.random.random((5,5))) ### ORIGINAL

		#### EXPERIMENTAL MOVE 
		im.set_data(flux[0]) ### ORIGINAL 
		#im.set_data(flux[0] / np.nanmedian(flux, axis=0))
		#im.set_data(np.nanmedian(flux, axis=0))  #### DOESN'T SEEM TO DO ANYTHING DIFFERENT. 
		fig.colorbar(im)


		ax.invert_xaxis() #### because RA increases to the left (EAST)	

		#ax.scatter(objx_pix, objy_pix, color='red', marker='*', s=100)
		ax.scatter(objx_pix, objy_pix, color='red', marker='*', s=100, label='obj')
		ax.scatter(refx_pix, refy_pix, color='green', marker='s', s=100, label='ref')

		try:
			ax.set_xticks(ticks=xticks, labels=xtick_coords)
			ax.set_yticks(ticks=yticks, labels=ytick_coords)
		except:
			ax.set_xticks(ticks=xticks)
			ax.set_yticks(ticks=yticks)
			ax.set_xticklabels(labels=xtick_coords)
			ax.set_yticklabels(labels=ytick_coords)		


		ax.set_xlabel('RA [deg]')
		ax.set_ylabel('Dec [deg]')


		return [im]

	# animation function.  This is called sequentially
	def animate(i):

		#a=im.get_array() #### ORIGINAL
		#a=a*np.exp(-0.001*i)    # exponential decay of the values ### ORIGINAL 
		next_frame = flux[i]
		#next_frame = flux[i] / np.nanmedian(flux, axis=0)

		#im.set_array(a) #### ORIGINAL
		ax.set_title('BKJD = '+str(round(times[i],2)))
		#ax.scatter(refx_pix, refy_pix, color='red', marker='*', s=100)
		#ax.scatter(objx_pix, objy_pix, color='red', marker='*', s=100, label='obj')
		#ax.scatter(refx_pix, refy_pix, color='green', marker='s', s=100, label='ref')
		"""
		try:
			ax.set_xticks(ticks=xticks, labels=xtick_coords)
			ax.set_yticks(ticks=yticks, labels=ytick_coords)
		except:
			ax.set_xticks(ticks=xticks)
			ax.set_yticks(ticks=yticks)
			ax.set_xticklabels(labels=xtick_coords)
			ax.set_yticklabels(labels=ytick_coords)		


		ax.set_xlabel('RA [deg]')
		ax.set_ylabel('Dec [deg]')
		"""
		im.set_array(next_frame)
		return [im]


	# call the animator.  blit=True means only re-draw the parts that have changed.
	#anim = animation.FuncAnimation(fig=fig, func=animate, fargs=(times, flux), init_func=init,
	#                               frames=flux.shape[0], interval=20, blit=True)

	anim = animation.FuncAnimation(fig=fig, func=animate, init_func=init, frames=flux.shape[0], repeat=False, interval=20, blit=True)
	#anim = animation.FuncAnimation(fig=fig, func=animate, frames=flux.shape[0], repeat=False, interval=20, blit=True)	

	# save the animation as an mp4.  This requires ffmpeg or mencoder to be
	# installed.  The extra_args ensure that the x264 codec is used, so that
	# the video can be embedded in html5.  You may need to adjust this for
	# your system: for more information, see
	# http://matplotlib.sourceforge.net/api/animation_api.html

	animation_filename = filename[:10]+'_FFI_animation.mp4'
	if save_animation == True:
		#anim.save(animation_filename, fps=30, extra_args=['-vcodec', 'libx264'])
		anim.save(animation_filename, fps=30)

	plt.show()


	if return_arrays == True:
		return times, flux 







