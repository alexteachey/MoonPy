from __future__ import division
import os
import numpy as np
import matplotlib.pyplot as plt 
from astropy.io import fits 



#### this will open the PRF file for convolution with a point source to overplot the appropriate aperture 
moonpydir = os.path.realpath(__file__)
moonpydir = moonpydir[:moonpydir.find('/mp_tess_prf.py')]
tess_prf_dir = moonpydir+'/TESS_PRF'



def get_PRF(camera, ccd, target_row, target_col, sector, show_plots=False):
	#### PRF files are sorted by camera number, ccd number, row number (y-value) and column number (x-value).
	##### NOTE THAT THERE IS NOT A PIXEL FOR EVERY SINGLE COORDINATE -- SO WE NEED TO INTERPOLATE BETWEEN TWO OF THEM
	##### ALSO NOTE THAT THERE ARE 9 PRF PIXELS PER TESS PIXEL.
	##### there are different PRFs for sectors 1-3 and sectors 4+

	#### find the appropriate file:
	if sector <= 3:
		sector_dir = tess_prf_dir+'/start_s0001'
	elif sector >= 4:
		sector_dir = tess_prf_dir+'/start_s0004'

	#### now establish the camera directory
	camera_dir = sector_dir+'/tess_prf_camera_'+str(camera)
	print('camera_dir: ', camera_dir)

	#### now establish the ccd directory
	ccd_dir = camera_dir+'/cam'+str(camera)+'_ccd'+str(ccd)
	print('ccd_dir: ', ccd_dir)

	prf_files = os.listdir(ccd_dir)
	rows, columns = [], []

	file_prefix = prf_files[0][:prf_files[0].find('row')]
	print('file_prefix = ', file_prefix)
	print(' ')

	for prf_file in prf_files:
		prf_row = prf_file[prf_file.find('row')+3:prf_file.find('row')+7] #### string! preserves leading zeros
		prf_col = prf_file[prf_file.find('col')+3:prf_file.find('col')+7] #### string! preserves leading zeros

		rows.append(prf_row) 
		columns.append(prf_col)

	rows, columns = np.array(rows), np.array(columns)
	unique_rows, unique_columns = np.sort(np.unique(rows)), np.sort(np.unique(columns))

	#### find the two closest rows to target row, and two closest columns to target column. Those will be the four corners.
	closest_rows = np.sort(unique_rows[np.argsort(np.abs(target_row - unique_rows.astype(int)))[:2]]) #### sorting still works even though they're strings! 
	closest_columns = np.sort(unique_columns[np.argsort(np.abs(target_col - unique_columns.astype(int)))[:2]]) ### sorting still works even though they're strings!

	closest_rows_int = closest_rows.astype(int)
	closest_columns_int = closest_columns.astype(int) 

	print('closest_rows_int: ', closest_rows_int)
	print('closest_columns_int: ', closest_columns_int)

	print('target_row = ', target_row)
	print('closest_rows = ', closest_rows)
	print(' ')
	print('target_col = ', target_col)
	print('closest_columns = ', closest_columns)
	print(' ')


	#### each row will be used twice, and each column will be used twice
	four_corner_files = []
	four_corner_prfs = []
	four_corner_uncs = []

	for ncr, cr in enumerate(closest_rows):
		for ncc, cc in enumerate(closest_columns):

			four_corner_filename = file_prefix+'row'+str(cr)+'-col'+str(cc)+'.fits'
			four_corner_files.append(four_corner_filename)


			#### open it up
			fcf_fits = fits.open(ccd_dir+'/'+four_corner_filename)
			fcf_prf = fcf_fits[0].data
			fcf_prf_unc = fcf_fits[1].data

			four_corner_prfs.append(fcf_prf)
			four_corner_uncs.append(fcf_prf_unc)

			### create the interpolation grids
			if (ncr == 0) and (ncc == 0):
				#### just do it once!

				dim1, dim2 = fcf_prf.shape[0], fcf_prf.shape[1]
				first_row_interp_grid = np.zeros(shape=(dim1, dim2))
				second_row_interp_grid = np.zeros(shape=(dim1, dim2))
				final_interp_grid = np.zeros(shape=(dim1, dim2))

			if (show_plots == True) or (show_plots == 'y'):
				plt.imshow(fcf_prf, origin='lower', interpolation='none')
				plt.title(four_corner_filename)
				plt.show()


	four_corner_files = np.array(four_corner_files)
	four_corner_prfs = np.array(four_corner_prfs)
	four_corner_uncs = np.array(four_corner_uncs)

	print("four_corner_prfs[0].shape = ", four_corner_prfs[0].shape)


	#### to interpolate -- we want to first interpolate between the lower row based on the target column
	#### then interpolate the upper row based on the target column,
	#### the interpolate between those two resulting images based on the target row 

	### because of the loop above, the first two four_corner_prfs are in the first row, 
	#### the second two four_corner_prfs are in the second row
	first_row_interp_row = int(closest_rows[0])
	second_row_interp_row = int(closest_rows[1])



	for i in np.arange(0,dim1,1):
		for j in np.arange(0,dim2,1):
			### np.interp takes (xvalue_at_interp_point, xvals, and yvals)
			first_row_interp_grid[i][j] = np.interp(x=int(target_col), xp=(closest_columns_int[0], closest_columns_int[1]), fp=(four_corner_prfs[0][i][j], four_corner_prfs[1][i][j]) )

			#### do the same thing with the second row
			second_row_interp_grid[i][j] = np.interp(x=int(target_col), xp=(closest_columns_int[0], closest_columns_int[1]), fp=(four_corner_prfs[2][i][j], four_corner_prfs[3][i][j]) )


			#### now interp between these two!
			final_interp_grid[i][j] = np.interp(x=int(target_row), xp=(first_row_interp_row, second_row_interp_row), fp=(first_row_interp_grid[i][j], second_row_interp_grid[i][j]) )


	##### see if it worked!
	if (show_plots == True) or (show_plots == 'y')
		plt.imshow(first_row_interp_grid, origin='lower', interpolation='none')
		plt.title('first row')
		plt.show()

		plt.imshow(second_row_interp_grid, origin='lower', interpolation='none')
		plt.title('second row')
		plt.show()

		plt.imshow(final_interp_grid, origin='lower', interpolation='none')
		plt.title('final interp grid')
		plt.show()








