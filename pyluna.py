from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os
from mp_tools import * 
import socket
from astropy.constants import G
import time


moonpydir = os.path.realpath(__file__)
moonpydir = moonpydir[:moonpydir.find('/pyluna.py')]

plt.rcParams["font.family"] = 'serif'

hostname = socket.gethostname()




### THIS CODE IS THE MASTER LUNA INTERFACE. ALL LUNA QUERIES THROUGH ME.
### probably should be written as a massive function that can be called.


### DIRECTORIES
#### NEW FUNCTIONALITY -- ALLOWS YOU TO RUN MULTIPLE LUNA runs at once -- your LUNAdir will be based in your CURRENT WORKING DIRECTORY!!!!
master_LUNAdir = moonpydir+'/LUNA'
local_LUNAdir = os.getcwd()+'/LUNA'
if (master_LUNAdir != local_LUNAdir):
	if os.path.exists(local_LUNAdir) == False:
		print('COPYING '+master_LUNAdir+' to '+local_LUNAdir)
		#### COPY master_LUNAdir locally -- that way you can make global changes within the master without doing ad-hoc changes.
		os.system('cp -r '+master_LUNAdir+' '+local_LUNAdir)
	
	LUNAdir = local_LUNAdir
	outputdir = LUNAdir+'/output'
	if os.path.exists(outputdir) == False:
		os.system('mkdir '+LUNAdir+'/output')


"""
#### make the make.sh file DIRECTORY SPECIFIC!
make_orig = open(LUNAdir+'/make_ORIG.sh', mode='r')
make_copy = open(LUNAdir+'/make_copy.sh', mode='w')

for nline,line in enumerate(make_orig):
	linesplit = line.split()
	if linesplit[0] == "$compiler":
		last_entry = LUNAdir+'/'
"""



"""
HERE IS THE MASTER FUNCTION! THIS IS WHAT YOU WILL CALL TO GENERATE THE MOONS YOU WANT.

"""

#def run_LUNA(tau0, nepochs, time_from_midtime_days, cadence_minutes, noise_ppm, star_params, plan_params, sat_params, munit='kg', runit='meters', ang_unit='radians', add_noise='y', show_plots='n', print_params='n', binned_output='y'):
#	assert len(star_params) == 4
#	assert len(plan_params) == 4
#	assert len(sat_params) == 6

def prepare_files(all_times, ntaus, nparam, nparamorig):
	#### ntaus is the number of transit times you will fit!
	#### if you are fitting TTVs, ntaus = the number of transits.
	#### if you are not fitting TTVs, ntaus = 1. 

	seriesP_file = open(LUNAdir+'/seriesP.jam', mode='w')

	for at in all_times:
		seriesP_file.write(str(at)+'\t1.\t0.001\t1\n')
	seriesP_file.close()


	### need to update plotit.f90 to make nplen = len(all_times)
	#plotitf90 = open(LUNAdir+'/plotit.f90', mode='r')
	#plotitf90_update = open(LUNAdir+'/plotit_update.f90', mode='w')
	plotitf90 = open(LUNAdir+'/plotit.f90', mode='r')
	plotitf90_update = open(LUNAdir+'/plotit_update.f90', mode='w')	

	for nline, line in enumerate(plotitf90):
		try:
			if line.split()[0] == 'nplen':
				newline = ' nplen = '+str(len(all_times))+'\n'

			elif line.split()[0] == 'nparam':
				newline = ' nparam = '+str(nparam)+'\n'
			elif line.split()[0] == 'OOTlen':
				newline = ' OOTlen = '+str(ntaus)+'\n'
			elif line.split()[0] == 'nparamorig':
				newline = ' nparamorig = '+str(nparamorig)+'\n'
			elif line.split()[0] == 'taulen':
				newline = ' taulen = '+str(ntaus)+'\n'
			else:
				newline = line 
		except:
			newline = line
		plotitf90_update.write(newline)
	plotitf90.close()
	plotitf90_update.close()
	#os.system('mv '+LUNAdir+'/plotit_update.f90 '+LUNAdir+'/plotit.f90')
	os.system('mv '+LUNAdir+'/plotit_update.f90 '+LUNAdir+'/plotit.f90')	

	print('sh '+LUNAdir+'/make.sh')
	#### NEED TO CHANGE DIRECTORY TO LUNAdir to make this file and have everything seen. Then change back.
	os.system('cd '+LUNAdir+' ; sh '+LUNAdir+'/make.sh')

	#os.system('sh make.sh')	





#def run_LUNA(all_times, tau0, Rstar, Mstar q1, q2, RpRstar, bplan, Pplan, RsatRp, MsatMp, sat_sma, sat_inc, sat_phase, sat_omega, cadence_minutes=29.42, noise_ppm=None, munit='kg', runit='meters', ang_unit='radians', add_noise='n', show_plots='n', print_params='n', binned_output='n'):
def run_LUNA(all_times, RpRstar, rhostar, bplan, Pplan, tau0, q1, q2, Psat=None, rhoplan=None, sat_sma=None, sat_phase=None, sat_inc=None, sat_omega=None, MsatMp=None, RsatRp=None, model="M", tau1=None, tau2=None, tau3=None, tau4=None, tau5=None, tau6=None, cadence_minutes=29.42, noise_ppm=None, munit='kg', runit='meters', ang_unit='radians', add_noise='n', show_plots='n', print_params='n', binned_output='n', **kwargs):
	#Rstar, Mstar, q1, q2 = star_params
	#Rplan, Mplan, bplan, Pplan = plan_params
	#Rsat, Msat, sat_sma, sat_inc, sat_phase, sat_omega = sat_params

	### have to make the file first! (Didn't use to be a problem).

	#### make sure units make sense
	#if runit == 'meters':
	#	assert (Rstar > 1e4) and (Rplan > 1e3) and (Rsat > 0) and (sat_sma < 1000)
	#else:
	#	pass ### build something ehre
	#if munit == 'kg':
	#	assert (Mstar > 1e27) and (Mplan > 1e20) and (Msat > 0)
	#lse:
	#	pass

	if ang_unit=='radians':
		#assert (sat_inc <= 2*np.pi) and (sat_phase <= 2*np.pi) and (sat_omega <= 2*np.pi)
		pass
	elif ang_unit=='degrees':
		#assert (sat_inc <= 360) and (sat_phase <= 360) and (sat_omega <= 360)
		### convert into radians!
		sat_inc = deg2rad(sat_inc)
		sat_phase = deg2rad(sat_phase)
		sat_omega = deg2rad(sat_omega)



	if rhoplan == None:
		### calculate it! 
		rhoplan = mp_tools.density_from_orbit(sat_sma, Psat, in_unit='days', out_unit='mks')


	#assert (bplan >= 0) and (bplan <= 2)

	### make the necessary conversions
	#RpRstar = Rplan / Rstar 
	#rhostar = mp_tools.density_conversion(Mstar, Rstar)
	#rhoplan = mp_tools.density_conversion(Mplan, Rplan)
	#rhosat = mp_tools.density_conversion(Msat, Rsat)
	#MsatMp = Msat / Mplan 
	#RsatRp = Rsat / Rplan 

	### misc conversions
	native_cadence = 0.2942439984 ### minutes
	native_cadence_days = native_cadence / (60 * 24) ### native for seriesP.jam
	cadence_days = cadence_minutes / (60 * 24)
	#epoch_midtimes = np.linspace(tau0, tau0+((nepochs-1)*Pplan), nepochs)

	#print('generating the input file...')
	### now you have to generate the input file!
	#input_file = open(LUNAdir+'/inputs.jam', mode='w')
	input_file = open(LUNAdir+'/inputs.jam', mode='w')	
	### inputs are 1) Rp/Rstar, rhostar, impact, Pplan, tau0, q1, q2, rho_plan, asp, phi_s, is, )s, Msp, Rsp.
	if 'e' in str(RpRstar):
		input_file.write(str('%.7f' % RpRstar)+'D0\n') #Rp(1)
	else:
		input_file.write(str(round(RpRstar,7))+'D0\n') #Rp(1)

	if 'e' in str(rhostar):
		input_file.write(str('%.7f' % rhostar)+'D0\n') #Rp(2)
	else:
		input_file.write(str(round(rhostar,7))+'D0\n') #Rp(2)

	if 'e' in str(bplan):
		input_file.write(str('%.7f' % bplan)+'D0\n') #Rp(3)
	else:
		input_file.write(str(round(bplan,7))+'D0\n') #Rp(3)

	if 'e' in str(Pplan):
		input_file.write(str('%.7f' % Pplan)+'D0\n') #Rp(4)
	else:
		input_file.write(str(round(Pplan,7))+'D0\n') #Rp(4)

	if 'e' in str(tau0):
		input_file.write(str('%.7f' % tau0)+'D0\n') #Rp(5)
	else:
		input_file.write(str(round(tau0,7))+'D0\n') #Rp(5)

	if 'e' in str(q1):
		input_file.write(str('%.7f' % q1)+'D0\n') #Rp(6)
	else:
		input_file.write(str(round(q1,7))+'D0\n') #Rp(6)

	if 'e' in str(q2):
		input_file.write(str('%.7f' % q2)+'D0\n') #Rp(7)
	else:
		input_file.write(str(round(q2,7))+'D0\n') #Rp(7)

	#if (model == 'M') or (model == 'Z'):

	if 'e' in str(rhoplan):
		input_file.write(str('%.7f' % rhoplan)+'D0\n') #Rp(8)
	else:
		input_file.write(str(round(rhoplan,7))+'D0\n') #Rp(8)

	if 'e' in str(sat_sma):
		input_file.write(str('%.7f' % sat_sma)+'D0\n') #Rp(9)
	else:
		input_file.write(str(round(sat_sma,7))+'D0\n') #Rp(9)

	if 'e' in str(sat_phase):
		input_file.write(str('%.7f' % sat_phase)+'D0\n') #Rp(10)
	else:
		input_file.write(str(round(sat_phase,7))+'D0\n') #Rp(10)

	if 'e' in str(sat_inc):
		input_file.write(str('%.7f' % sat_inc)+'D0\n') #Rp(11)
	else:
		input_file.write(str(round(sat_inc,7))+'D0\n') #Rp(11)

	if 'e' in str(sat_omega):
		input_file.write(str('%.7f' % sat_omega)+'D0\n') #Rp(12)
	else:
		input_file.write(str(round(sat_omega,7))+'D0\n') #Rp(12)

	if 'e' in str(MsatMp):
		input_file.write(str('%.7f' % MsatMp)+'D0\n') # Rp(13)
	else:
		input_file.write(str(round(MsatMp,7))+'D0\n') # Rp(13)
	if 'e' in str(RsatRp):
		input_file.write(str('%.7f' % RsatRp)+'D0\n') # Rp(14)
	else:
		input_file.write(str(round(RsatRp,7))+'D0\n') # Rp(14)	


	#### FITTING INDIVIDUAL TAUS!
	if model == 'T':
		if tau1 != None:
			if 'e' in str(tau1):
				input_file.write(str('%.7f' % tau1)+'D0\n') # Rp(15)
			else:
				input_file.write(str(round(tau1,7))+'D0\n') # Rp(15)	
		#### FITTING INDIVIDUAL TAUS!
		if tau2 != None:
			if 'e' in str(tau2):
				input_file.write(str('%.7f' % tau2)+'D0\n') # Rp(16)
			else:
				input_file.write(str(round(tau2,7))+'D0\n') # Rp(16)	

		#### FITTING INDIVIDUAL TAUS!
		if tau3 != None:
			if 'e' in str(tau3):
				input_file.write(str('%.7f' % tau3)+'D0\n') # Rp(17)
			else:
				input_file.write(str(round(tau3,7))+'D0\n') # Rp(17)	

		#### FITTING INDIVIDUAL TAUS!
		if tau4 != None:
			if 'e' in str(tau4):
				input_file.write(str('%.7f' % tau4)+'D0\n') # Rp(18)
			else:
				input_file.write(str(round(tau4,7))+'D0\n') # Rp(18)

		#### FITTING INDIVIDUAL TAUS!
		if tau5 != None:
			if 'e' in str(tau5):
				input_file.write(str('%.7f' % tau5)+'D0\n') # Rp(19)
			else:
				input_file.write(str(round(tau5,7))+'D0\n') # Rp(19)

		#### FITTING INDIVIDUAL TAUS!
		if tau6 != None:
			if 'e' in str(tau1):
				input_file.write(str('%.7f' % tau6)+'D0\n') # Rp(20)
			else:
				input_file.write(str(round(tau6,7))+'D0\n') # Rp(20)			


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
		if (model == 'M') or (model == "Z"):
			print("planet density [kg / m^3] = ", rhoplan)
			print("sat_sma = [Rp] ", sat_sma)
			print("sat_phase = ", sat_phase)
			print("sat_inc = ", sat_inc)
			print("sat_omega = ", sat_omega)
			print("Msat / Mp = ", MsatMp)
			print("Rsat / Rp = ", RsatRp)
		print(" ")

	### now it's time to run plotit!
	print("calling plotit.")
	os.system('cd '+LUNAdir+' ; ./plotit')


	### now you should have light curves... load them and plot them
	#moon_file = np.genfromtxt(outputdir+'/PRI_full.0.jam')
	moon_file = np.genfromtxt(LUNAdir+'/output/PRI_full.0.jam')
	output_times, output_fluxes = moon_file.T[0], moon_file.T[1]

	#if np.all(output_fluxes == 1):
	#	raise Exception('flat light curve.')

	#if len(times) != len(output_times):
	#	raise Exception('len(data) != len(model).')

	#print("pyluna times = ", output_times)
	#print('pyluna fluxes = ', output_fluxes)


	"""
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
		lc_binidxs = np.digitize(output_times, lc_bins)
		for nlcbin, lcbin in enumerate(lc_bins):
			### grab the fluxes for this bin
			bin_fluxes = output_fluxes[np.where(lc_binidxs == nlcbin)]
			avg_bin_fluxes = np.nanmean(bin_fluxes)
			binned_fluxes.append(avg_bin_fluxes)
		binned_fluxes = np.array(binned_fluxes)



	if show_plots == 'y':
		if binned_output == 'y':
			plt.scatter(lc_bins, binned_fluxes, facecolors='LightCoral', edgecolors='k', s=20)
		elif binned_output == 'n':
			plt.scatter(output_times, output_fluxes, facecolors='LightCoral', edgecolors='k', s=20)

		plt.xlabel('Time')
		plt.ylabel('Flux')
		plt.show()

		if add_noise == 'y':
			if binned_output == 'y':
				noisy_fluxes = np.random.normal(loc=output_fluxes, scale=noise_ppm*1e-6)
			elif binned_output == 'n':
				noisy_fluxes = np.random.normal(loc=output_fluxes, scale=noise_ppm*1e-6)
			
			plt.scatter(output_times, noisy_fluxes, facecolors='LightCoral', edgecolors='k', s=20)
			plt.xlabel('Time')
			plt.ylabel('Flux')
			plt.show()
	"""

	"""
	if add_noise == 'n':
		if binned_output == 'y':
			return lc_bins, binned_fluxes
		elif binned_output == 'n':
			return output_times, output_fluxes

	elif add_noise == 'y':
		if binned_output == 'y':
			return lc_bins, noisy_fluxes
		elif binned_output == 'n':
			return output_times, noisy_fluxes 
	"""

	return output_times, output_fluxes 


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