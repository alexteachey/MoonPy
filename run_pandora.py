from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os
#from mp_tools import * 
import mp_tools
import socket
from astropy.constants import G
import time
try:
	import pandoramoon as pandora
	from pandoramoon.helpers import ld_convert, ld_invert 
except:
	print("could not import pandora. You ca 'pip install pandoramoon' to rectify this. ")


moonpydir = os.path.realpath(__file__)
moonpydir = moonpydir[:moonpydir.find('/run_pandora.py')]

plt.rcParams["font.family"] = 'serif'

hostname = socket.gethostname()


### THIS CODE IS THE MASTER PANDORA INTERFACE. ALL PANDORA QUERIES THROUGH ME.
### probably should be written as a massive function that can be called.


#def run_Pandora(all_times, Rstar_meters, Pplan, apRstar, RpRstar, bplan, tau0, Mplan_kg, eccplan, Tdur_days, RsRstar, Psat, sat_phase, sat_inc, sat_omega, Msat_kg, q1, q2, input_ang_unit='degrees', cadence_minutes=29.42, t0_offset=0.01, sat_Omega=0, print_params=True, **kwargs):
def run_Pandora(all_times, nepochs, R_star, per_bary, a_bary, r_planet, b_bary, t0_bary_offset, M_planet, r_moon, per_moon, tau, Omega_moon, i_moon, M_moon, q1, q2, t0_bary, Tdur_days, ecc_bary=0, w_moon=0, input_ang_unit='degrees', cadence_minutes=29.42, t0_offset=0.01, sat_Omega=0, print_params=True, **kwargs):

	#### removed t0_bary 
	#def run_Pandora(all_times, R_star, per_bary, a_bary, r_planet, b_bary, t0_bary_offset, M_planet, r_moon, per_moon, tau, Omega_moon, i_moon, M_moon, q1, q2, Tdur_days, nepochs, ecc_bary=0, w_moon=0, input_ang_unit='degrees', cadence_minutes=29.42, t0_offset=0.01, sat_Omega=0, print_params=True, **kwargs):

	if input_ang_unit=='radians':
		#### PANDORA WANTS DEGREES
		sat_inc = rad2deg(sat_inc)
		sat_phase = rad2deg(sat_phase)
		sat_omega = rad2deg(sat_omega)


	elif input_ang_unit=='degrees':
		pass


	### misc conversions
	cadence_days = cadence_minutes / (60 * 24)



	#### 15 parameters in all -- 
	#### see here: http://localhost:8888/notebooks/injection_retrieval_simple_ultranest.ipynb#:~:text=Create%20planet%2Bmoon%20model 

	params = pandora.model_params()

	#### VARIABLE PARAMETERS (FOR FITTING)
	u1, u2 = ld_convert(q1, q2) 
	
	params.R_star = float(R_star) #### FIT PARAM #0 
	params.per_bary = float(per_bary) #### PARAM #1 Pplan [days]
	params.a_bary = float(a_bary)  #### PARAM #2 
	params.r_planet = float(r_planet) #### PARAM #3 
	params.b_bary = float(b_bary) #### PARAM # 4
	params.t0_bary_offset = float(t0_bary_offset) #### PARAM #5 what is this? 
	params.M_planet = float(M_planet) #### PARAM #6 [kg]
	#### moon parameters
	params.r_moon = float(r_moon) #### PARAM #7 -- satellite radius divided by stellar radius
	params.per_moon = float(per_moon) #### PARAM #8 -- need to define above
	params.tau_moon = float(tau) #### PARAM #9-- must be between zero and one -- I think this is the phase...
	params.Omega_moon = float(Omega_moon) #### PARAM # 10 -- longitude of the ascending node??? between 0 and 180 degrees 
	params.i_moon = float(i_moon) #### PARAM #11 -- between 0 and 180 degrees
	#params.M_moon = MsatMp * Mplan_kg #### PARAM # 12 need to define Mp! [kg]
	params.M_moon = float(M_moon) 
	params.u1 = float(u1) #### PARAM #13 need to define above!
	params.u2 = float(u2) #### PARAM #14 need to define above!


	#### FIXED PARAMETERS -- BUT I"m NOT SURE WHY ... 
	params.t0_bary = float(t0_bary) #### FIXED?
	#params.ecc_bary = float(ecc_bary) #### FIXED? -- need to define above!	
	#params.w_moon = float(w_moon) #### degrees

	#### other inputs
	params.epochs = nepochs #### needs to be defined above!
	params.epoch_duration = 20*Tdur_days #### need to be defined above
	params.cadences_per_day = 48 
	#params.epoch_distance = Pplan 
	params.epoch_distance = per_bary 
	params.supersampling_factor = 1
	params.occult_small_threshold = 0.1 ### between 0 and 1 -- what is this?
	params.hill_sphere_threshold = 1.2 #### what does this mean?


	pdtime = pandora.time(params).grid()
	pdmodel = pandora.moon_model(params)

	#total_flux, planet_flux, moon_flux = pdmodel.light_curve(pdtime)
	#total_flux, planet_flux, moon_flux = pdmodel.light_curve(all_times)
	total_flux, planet_flux, moon_flux = pdmodel.light_curve(all_times)


	if print_params == 'y':

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


	#return output_times, output_fluxes 
	return total_flux, planet_flux, moon_flux







