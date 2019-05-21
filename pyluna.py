from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os
import mp_tools 


### THIS CODE IS THE MASTER LUNA INTERFACE. ALL LUNA QUERIES THROUGH ME.
### probably should be written as a massive function that can be called.


### DIRECTORIES
LUNAdir = '/Users/hal9000/Documents/Software/MoonPy/LUNA'
outputdir = LUNAdir+'/output'


"""
HERE IS THE MASTER FUNCTION! THIS IS WHAT YOU WILL CALL TO GENERATE THE MOONS YOU WANT.

"""

#def run_LUNA(tau0, nepochs, time_from_midtime_days, cadence_minutes, noise_ppm, star_params, plan_params, sat_params, munit='kg', runit='meters', ang_unit='radians', add_noise='y', show_plots='n', print_params='n', binned_output='y'):
#	assert len(star_params) == 4
#	assert len(plan_params) == 4
#	assert len(sat_params) == 6


def run_LUNA(times, tau0, Rstar, Mstar, q1, q2, Rplan, Mplan, bplan, Pplan, Rsat, Msat, sat_sma, sat_inc, sat_phase, sat_omega, cadence_minutes=29.42, noise_ppm=None, munit='kg', runit='meters', ang_unit='radians', add_noise='y', show_plots='n', print_params='n', binned_output='n'):
	#Rstar, Mstar, q1, q2 = star_params
	#Rplan, Mplan, bplan, Pplan = plan_params
	#Rsat, Msat, sat_sma, sat_inc, sat_phase, sat_omega = sat_params

	### make sure units make sense
	if runit == 'meters':
		assert (Rstar > 1e4) and (Rplan > 1e3) and (Rsat > 0) and (sat_sma < 1000)
	else:
		pass ### build something ehre
	if munit == 'kg':
		assert (Mstar > 1e27) and (Mplan > 1e20) and (Msat > 0)
	else:
		pass

	if ang_unit=='radians':
		assert (sat_inc <= 2*np.pi) and (sat_phase <= 2*np.pi) and (sat_omega <= 2*np.pi)
	elif ang_unit=='degrees':
		assert (sat_inc <= 360) and (sat_phase <= 360) and (sat_omega <= 360)
		### convert into radians!
		sat_inc = deg2rad(sat_inc)
		sat_phase = deg2rad(sat_phase)
		sat_omega = deg2rad(sat_omega)

	assert (bplan >= 0) and (bplan <= 2)

	### make the necessary conversions
	RpRstar = Rplan / Rstar 
	rhostar = density_conversion(Mstar, Rstar)
	rhoplan = density_conversion(Mplan, Rplan)
	rhosat = density_conversion(Msat, Rsat)
	MsatMp = Msat / Mplan 
	RsatRp = Rsat / Rplan 

	### misc conversions
	native_cadence = 0.2942439984 ### minutes
	native_cadence_days = native_cadence / (60 * 24) ### native for seriesP.jam
	cadence_days = cadence_minutes / (60 * 24)
	#epoch_midtimes = np.linspace(tau0, tau0+((nepochs-1)*Pplan), nepochs)


	### now the first thing you need to do is generate the seriesP.jam file... this is read in to create the light curves.
	seriesP_file = open(LUNAdir+'/seriesP.jam', mode='w')
	#for nep, epoch_midtime in enumerate(epoch_midtimes):
	#	epoch_times = np.arange(epoch_midtime-time_from_midtime_days, epoch_midtime+time_from_midtime_days+native_cadence_days, native_cadence_days)
	#	if nep == 0:
	#		all_times = epoch_times
	#	else:
	#		all_times = np.concatenate((all_times, epoch_times))
	all_times = times 
	#	for transit_time in epoch_times:
	#		seriesP_file.write(str(transit_time)+'\t1.\t0.001\t1\n')
	for at in all_times:
		seriesP_file.write(str(at)+'\t1.\t0.001\t1\n')
	seriesP_file.close()


	### need to update plotit.f90 to make nplen = len(all_times)
	plotitf90 = open(LUNAdir+'/plotit.f90', mode='r')
	plotitf90_update = open(LUNAdir+'/plotit_update.f90', mode='w')

	for nline, line in enumerate(plotitf90):
		try:
			if line.split()[0] == 'nplen':
				newline = ' nplen = '+str(len(all_times))+'\n'
			else:
				newline = line 
		except:
			newline = line
		plotitf90_update.write(newline)
	plotitf90.close()
	plotitf90_update.close()
	os.system('mv '+LUNAdir+'/plotit_update.f90 '+str(LUNAdir)+'/plotit.f90')



	### now you have to generate the input file!
	input_file = open(LUNAdir+'/inputs.jam', mode='w')
	### inputs are 1) Rp/Rstar, rhostar, impact, Pplan, tau0, q1, q2, rho_plan, asp, phi_s, is, )s, Msp, Rsp.
	input_file.write(str(round(RpRstar,7))+'D0\n') #Rp(1)
	input_file.write(str(round(rhostar,7))+'D0\n') #Rp(2)
	input_file.write(str(round(bplan,7))+'D0\n') #Rp(3)
	input_file.write(str(round(Pplan,7))+'D0\n') #Rp(4)
	input_file.write(str(round(tau0,7))+'D0\n') # Rp(5)
	input_file.write(str(round(q1,7))+'D0\n') #Rp(6) ### UPDATED FEBRUARY 18th, 2019!
	input_file.write(str(round(q2,7))+'D0\n') #Rp(7) ### UPDATED FEBRUARY 18th, 2019!
	input_file.write(str(round(rhoplan,7))+'D0\n') #Rp(8)
	input_file.write(str(round(sat_sma,7))+'D0\n') #Rp(9)
	input_file.write(str(round(sat_phase,7))+'D0\n') #Rp(10)
	input_file.write(str(round(sat_inc,7))+'D0\n') # Rp(11)
	input_file.write(str(round(sat_omega,7))+'D0\n') # RP(12)
	if 'e' in str(MsatMp):
		input_file.write(str('%.7f' % MsatMp)+'D0\n')
	else:
		input_file.write(str(round(MsatMp,7))+'D0\n') # Rp(13)
	if 'e' in str(RsatRp):
		input_file.write(str('%.7f' % RsatRp)+'D0\n')
	else:
		input_file.write(str(round(RsatRp,7))+'D0\n') # Rp(14)	

	input_file.close()


	if print_params == 'y':
		#print("nepochs = ", nepochs)
		#print("midtimes = ", epoch_midtimes)
		#print("native_cadence = ", native_cadence)
		#print("native_cadence_days = ", native_cadence_days)
		#print("cadence_minutes = ", cadence_minutes)
		#print("cadence_days = ", cadence_days)
		print(" ")
		print("Rp/Rstar = ", RpRstar)
		print("transit depth [ppm] = ", RpRstar**2 * 1e6)
		print("stellar density [kg / m^3] = ", rhostar)
		print("impact = ", bplan)
		print("Period [days] = ", Pplan)
		print("tau_0 [day] = ", tau0)
		print("q1,q2 = ", q1, q2)
		print("planet density [kg / m^3] = ", rhoplan)
		print("sat_sma = [Rp] ", sat_sma)
		print("sat_phase = ", sat_phase)
		print("sat_inc = ", sat_inc)
		print("sat_omega = ", sat_omega)
		print("Msat / Mp = ", MsatMp)
		print("Rsat / Rp = ", RsatRp)
		print(" ")

	### now it's time to run plotit!
	os.system('./plotit')


	### now you should have light curves... load them and plot them
	moon_file = np.genfromtxt(outputdir+'/PRI_full.0.jam')
	times, fluxes = moon_file.T[0], moon_file.T[1]

	if binned_output == 'y':
		binned_fluxes = []
		### bin the results!
		lc_bins = np.arange(np.nanmin(all_times)-(0.5*cadence_days), np.nanmax(all_times)+(0.5*cadence_days)+cadence_days, cadence_days)
		print("len(lc_bins) (pre-reduce) = ", len(lc_bins))
		### remove out of bounds bins
		bad_bin_idxs = []
		for nlcbin, lcbin in enumerate(lc_bins):
			if np.all(np.abs(lcbin - epoch_midtimes) > time_from_midtime_days):
				bad_bin_idxs.append(nlcbin)
		lc_bins = np.delete(lc_bins, bad_bin_idxs)
		print("len(lc_bins) (post-reduce) = ", len(lc_bins))




		print("lc_bins = ", lc_bins)
		lc_binidxs = np.digitize(times, lc_bins)
		for nlcbin, lcbin in enumerate(lc_bins):
			### grab the fluxes for this bin
			bin_fluxes = fluxes[np.where(lc_binidxs == nlcbin)]
			avg_bin_fluxes = np.nanmean(bin_fluxes)
			binned_fluxes.append(avg_bin_fluxes)
		binned_fluxes = np.array(binned_fluxes)



	if show_plots == 'y':
		if binned_output == 'y':
			plt.scatter(lc_bins, binned_fluxes, facecolors='LightCoral', edgecolors='k', s=20)
		elif binned_output == 'n':
			plt.scatter(times, fluxes, facecolors='LightCoral', edgecolors='k', s=20)

		plt.xlabel('Time')
		plt.ylabel('Flux')
		plt.show()

		if add_noise == 'y':
			if binned_output == 'y':
				noisy_fluxes = np.random.normal(loc=fluxes, scale=noise_ppm*1e-6)
			elif binned_output == 'n':
				noisy_fluxes = np.random.normal(loc=fluxes, scale=noise_ppm*1e-6)
			
			plt.scatter(times, noisy_fluxes, facecolors='LightCoral', edgecolors='k', s=20)
			plt.xlabel('Time')
			plt.ylabel('Flux')
			plt.show()


	if add_noise == 'n':
		if binned_output == 'y':
			return lc_bins, binned_fluxes
		elif binned_output == 'n':
			return times, fluxes

	elif add_noise == 'y':
		if binned_output == 'y':
			return lc_bins, noisy_fluxes
		elif binned_output == 'n':
			return times, noisy_fluxes 


"""
### for testing purposes
test_tau0 = 1325
test_nepochs = 12
test_time_from_tmid = 2
test_cadence = 29.42
test_noise = 50
test_star_params = np.array([eq_RSun, MSun, 0.4329845645336493, 0.3765240747629319]) ### Rstar, Mstar, q1, q2
test_plan_params = np.array([eq_RJup, MJup, 0.1, 34.5558]) ### Rplan, Mplan, bplan, Pplan
test_sat_params = np.array([0.5*eq_REarth, 0.1*MEarth, 10, 0, 0, 0]) ### Rsat, Msat, sat_sma, sat_inc, sat_phase, sat_omega 


### test it!

LUNA_times, LUNA_fluxes = run_LUNA(test_tau0, test_nepochs, test_time_from_tmid, test_cadence, test_noise, test_star_params, test_plan_params, test_sat_params, add_noise='y', show_plots='y', print_params='y')
"""