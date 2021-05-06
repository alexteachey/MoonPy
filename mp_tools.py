from __future__ import division
import numpy as np 
from astropy.constants import G, c, M_earth, M_jup, M_sun, R_earth, R_jup, R_sun, au 
import os
import pandas


moonpydir = os.path.realpath(__file__)
moonpydir = moonpydir[:moonpydir.find('/mp_tools.py')]

DKS_LDCs = pandas.read_csv(moonpydir+'/DKSing_LDCs_Kepler_Table2_format_mod.csv')
DKS_TeffK = np.array(DKS_LDCs['Teff_K']).astype(float)
DKS_Logg = np.array(DKS_LDCs['Logg']).astype(float)
DKS_MH = np.array(DKS_LDCs['M/H']).astype(float)
DKS_quada = np.array(DKS_LDCs['quad_a']).astype(float)
DKS_quadb = np.array(DKS_LDCs['quad_b']).astype(float)


def DKS_best_LDCmatch(Teff=np.nan, Logg=np.nan, MH=np.nan):
	#### find the best match in the limb-darkening catalog
	Teff_diffsq = (DKS_TeffK - Teff)**2
	Logg_diffsq = (DKS_Logg - Logg)**2
	MH_diffsq = (DKS_MH - MH)**2
	absolute_distances = np.sqrt(np.nansum(np.array([Teff_diffsq, Logg_diffsq, MH_diffsq]), axis=0))

	best_match_idx = np.nanargmin(absolute_distances)
	best_match_quada, best_match_quadb = DKS_quada[best_match_idx], DKS_quadb[best_match_idx]
	return best_match_quada, best_match_quadb  




### FROM ASTROPY
eq_RSun = R_sun.value ### meters
eq_RJup = R_jup.value ### meters
eq_Rearth = R_earth.value ### meters
MEarth = M_earth.value ### kg
MJup = M_jup.value ### kg
MSun = M_sun.value ### kg 

### FROM WOLFRAM
eq_RNep = 2.4622e7 ### meters
eq_RMoon = 1.7374e6 ### meters 
rhoEarth = 5515 ### kg / m^3
rhoMoon = 3344 ### kg / m^3
rhoNep = 1638 ### kg / m^3
rhoJup = 1326.2 ### kg / m^3
rhoSun = 1408 ### kg / m^3
MMoon = 7.3459e22 ### kg
MNep = 1.0241e26 ### kg
sma_moon = 3.844e8 ### meters



#### FUNCTIONS 
def effective_radius(density, mass):
	### CALCULATE AN EFFECTIVE RADIUS OF THE MOON FROM A PROVIDED MASS AND DENSITY -- DIFFERENT FROM EQUATORIAL!
	return ((3/4) * (mass/(np.pi * density)))**(1/3)

def RHill(sma_plan, m_star, m_plan): ### will work as long as both masses are in the same units!
	return sma_plan * ((m_plan / (3*m_star)))**(1/3)

def Kep3_pfroma(sma, m1, m2, sma_unit = 'meters', output_format='days', val_only='y'):
	### CALCULATE THE PERIOD BASED ON THE SEMIMAJOR AXIS
	native_solution = np.sqrt((sma**3 * 4 * np.pi**2)/(G*(m1+m2)))
	if output_format == 'days':
		output = native_solution / (24 * 60 * 60)
	else:
		output = native_solution
	if val_only == 'y':
		output = output.value 
	return output 

def Kep3_afromp(period, m1, m2, val_only='y', unit='days'):
	if unit=='days':
		period = period * (24 * 60 * 60) ### convert to seconds!
	### CALCULATE THE SEMIMAJOR AXIS BASED ON THE PERIOD
	numerator = period**2 * G * (m1 + m2)
	denominator = 4*np.pi**2
	sma = (numerator/denominator)**(1/3)
	if val_only == 'y':
		sma = sma.value
	return sma 

def mass_from_density(density, radius):
	### density and radius should have matching units!
	object_volume = (4/3) * np.pi * radius**3
	object_mass = density * object_volume
	return object_mass 

def inc_from_impact(impact, rstar, sma, unit='radians'):
	### units must be consistent!
	inclination = np.arccos((impact * rstar)/sma)
	if unit == 'degrees':
		inclination = (inclination * (180 / np.pi))
	else:
		inclination = inclination
	return inclination 


def impact_from_inc(inclination, rstar, sma, unit='degrees'):
	### units must be consistent!
	if unit == 'degrees':
		inclination = inclination * (np.pi / 180)
	impact = (sma * np.cos(inclination)) / rstar
	return impact 

def Tdur(period, Rstar, Rplan, impact, sma):
	### CALCULATE THE DURATION OF A PLANETARY TRANSIT.
	first_term = period/np.pi
	second_term_numerator = np.sqrt((Rstar+Rplan)**2 - (impact*Rstar)**2)
	second_term_denominator = sma
	return first_term * np.arcsin(second_term_numerator/second_term_denominator)


def quadsum(values):
	### SUM VALUES IN QUADRATURE
	return np.sqrt(np.sum(values**2))

def deg2rad(degrees):
	return degrees * (np.pi/180)

def rad2deg(radians):
	return radians * (180/np.pi)


def u1u2_to_q1q2(u1, u2): 
	### this function converts the standard, quadratic limb darkening coefficents to 
	### David Kipping's q1 and q2, which are the standard inputs for LUNA. 
	### see arXiv 1308.0009 equations 17 and 18 for this reparameterization
	q1 = (u1 + u2)**2
	q2 = (0.5 * u1) / (u1 + u2)
	return q1, q2 

def q1q2_to_u1u2(q1, q2):
	u1 = 2*q2*np.sqrt(q1)
	u2 = -2*(q2 - 0.5) * np.sqrt(q1)
	return u1, u2


def Rp_timescale(times, impact, sma_plan, Rplan, Rstar, Pplan, Tmid):
	### note that all time and size units must be the same!
	### this is equation 8 in Teachey et al 2018
	first_term = impact**2
	second_term = (sma_plan/Rstar)**2 - impact**2
	third_term_arg = ((2*np.pi)/Pplan) * (times - Tmid)
	third_term = (np.sin(third_term_arg))**2
	fourth_term = 1
	bracket_term = np.sqrt(first_term + (second_term * third_term)) - fourth_term 
	p = Rplan / Rstar
	final_answer = bracket_term / p 
	return final_answer 


def Roche(Rsat, Mplan, Msat):
	return Rsat * ((2*Mplan) / Msat)**(1/3)


def transit_SNR(depth, error_per_obs, transit_duration_minutes):
	### SIMPLE TRANSIT SNR CALCULATOR -- DOESN'T ACCOUNT FOR TRANSIT MORPHOLOGY.
	return (depth/error_per_obs) * np.sqrt((transit_duration_minutes / 30))


def density_conversion(mass, radius, munit='kg', runit='meters'):
	rho_obj = mass / ((4/3) * np.pi * radius**3)
	return rho_obj 


def density_from_orbit(a_over_R, Porbit, in_unit='days', out_unit='mks'):
	#### this computes a density based on the semimajor axis of the orbiter and the orbital period.
	##### SEE EQUATION 3 HERE: https://arxiv.org/pdf/1710.07293.pdf
	numerator = 3*np.pi * a_over_R**3
	if in_unit == 'days':
		Porbit_seconds = Porbit*24*60*60
	elif in_unit == 'hours':
		Porbit_seconds = Porbit*60*60
	denominator = G.value * Porbit_seconds**2
	density_kgm3 = numerator / denominator
	if out_unit == 'mks':
		return density_kgm3
	elif out_unit == 'cgs':
		density_gcm3 = density_kgm3 * 1000 * (1/100)**3
		return density_gcm3



def transit_SNR_integrator(times, fluxes, errors):
	### THIS CAN ONLY BE USED ON A FLAT MODEL! DOESN'T WORK FOR FLUXES!
	### make short times and fluxes are sorted
	timesort = np.argsort(times)
	times = times[timesort]
	fluxes = fluxes[timesort]

	### identify the first and last instances where flux < 1
	transit_idxs = np.where((fluxes < 1))[0]
	min_transit_idxs = np.where((fluxes == np.nanmin(fluxes)))[0]

	Tobs = np.nanmax(times) - np.nanmin(times)
	T1_idx, T2_idx, T3_idx, T4_idx = transit_idxs[0], min_transit_idxs[0], min_transit_idxs[-1], transit_idxs[-1]
	T1, T2, T3, T4 = times[T1_idx], times[T2_idx], times[T3_idx], times[T4_idx]
	T1_flux, T2_flux, T3_flux, T4_flux = fluxes[T1_idx], fluxes[T2_idx], fluxes[T3_idx], fluxes[T4_idx]
	T14, T23 = T4-T1, T3-T2

	### taken from Kipping and Sandford https://arxiv.org/pdf/1603.05662.pdf
	W = (T14+T23)/2 

	### in the limit where the depth << 1 and Tobs >> T_transit, we can right

	transit_dt = times[T1_idx+1] - times[T1_idx] ### assumes uniform cadence. ### UNITS OF DAYS!!!
	nobs_per_LC = cadence_days / transit_dt ### number of hires data points in the long cadence.
	depths = 1 - fluxes[transit_idxs]
	transit_depth = np.nanmax(depths)

	hires_precision = precision * np.sqrt(nobs_per_LC)
	single_point_SNRs = depths / hires_precision

	quadsum_single_point_SNRs = quadsum(single_point_SNRs)
	final_SNR = quadsum_single_point_SNRs 


	area_of_the_transit = np.sum(depths * transit_dt)

	SNR_sanity_check = (transit_depth / precision) * np.sqrt(area_of_the_transit / (transit_depth * cadence_days))

	if run_diagnostics == 'y':
		print("----------------")
		print("SNR params: ")
		print("Baseline [days] = ", Tobs)
		print("# obs in transit = ", T14/cadence_days)
		print("Duration (T14) [days] = ", T14)
		print("transit_depth [ppm] = ", transit_depth*1e6)
		print("stellar precision [ppm] = ", (precision*1e6))
		print("quadsum_single_point_SNRs = ", quadsum_single_point_SNRs)
		print("final_SNR = ", final_SNR)
		print("SNR_sanity_check = ", SNR_sanity_check)
		print("---------------")
		print(" ")

	return final_SNR 



def lc_fold(times, fluxes, errors, tau0, period, phase_offset=0.0):
	### this method will phase fold your light curve. 
	### first tau in the time series:
	first_tau = tau0

	fold_times = ((((np.hstack(times) - first_tau - 0.5*period - phase_offset*period) % period) / period)) ### yields the remainder!
	fold_times = fold_times - 0.5

	fold_fluxes = np.hstack(fluxes)
	fold_errors = np.hstack(errors)

	##### sort them
	fold_sort_idxs = np.argsort(fold_times)
	fold_times, fold_fluxes, folded_errors = fold_times[fold_sort_idxs], fold_fluxes[fold_sort_idxs], fold_errors[fold_sort_idxs]

	return fold_times, fold_fluxes, fold_errors


def nospaces(string):
	newstring = string.replace(' ', '')
	return newstring 


def DWstat(data, model):
	residuals = data - model
	residual_difs = [] 
	for i in np.arange(1,len(residuals),1):
		residual_dif = residuals[i] - residuals[i-1]
		residual_difs.append(residual_dif)
	residual_difs = np.array(residual_difs)
	DWnum = np.nansum(residual_difs**2)
	DWdenom = np.nansum(residuals**2)
	DW = DWnum / DWdenom
	return DW 



