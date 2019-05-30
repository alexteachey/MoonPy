from __future__ import division
import numpy as np
import matplotlib.pyplot as plt 
from mp_tools import mass_from_density, q1q2_to_u1u2, inc_from_impact, Kep3_afromp



def run_batman(all_times, RpRstar, Rstar, bplan, Pplan, tau0, q1, q2, long_peri=0, ecc=0, Mstar=None, Mplan=None, rhostar=None, rhoplan=None, cadence_minutes=29.42, noise_ppm=None, munit='kg', runit='meters', ang_unit='radians', add_noise='n', show_plots='n', print_params='n', binned_output='n', **kwargs):
	import batman
	#### initial calculations
	### you may supply planet masses OR densities!
	all_times = np.hstack(all_times)
	if Mstar == None:
		### you must have rhostar!
		Mstar = mass_from_density(rhostar, Rstar)
	if Mplan == None:
		Mplan = mass_from_density(rhoplan, RpRstar*Rstar)
	planet_sma = Kep3_afromp(Pplan, Mstar, Mplan)
	planet_sma_Rstar = planet_sma/Rstar 
	#print('planet_sma = ', planet_sma)
	#print('impact parameter = ', bplan)
	#print('inclination = ', inc_from_impact(bplan, Rstar, planet_sma, unit='degrees'))

	batman_params = batman.TransitParams()

	batman_params.t0 = tau0
	batman_params.per = Pplan ### in days
	batman_params.rp = RpRstar ### natively in stellar units
	batman_params.a = planet_sma_Rstar
	batman_params.inc = inc_from_impact(bplan, Rstar, planet_sma, unit='degrees')
	batman_params.ecc = ecc 
	batman_params.w = long_peri
	u1, u2 = q1q2_to_u1u2(q1, q2)
	batman_params.u = [u1, u2]
	batman_params.limb_dark = 'quadratic'

	batman_model = batman.TransitModel(batman_params, all_times)
	batman_fluxes = batman_model.light_curve(batman_params)

	if show_plots=='y':
		plt.plot(all_times, batman_fluxes)
		plt.show()

	return all_times, batman_fluxes 
