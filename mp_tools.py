from __future__ import division
import numpy as np 
#import mr_forecast ### needs to be somewhere readable!
from astropy.constants import G, c, M_earth, M_jup, M_sun, R_earth, R_jup, R_sun, au 

### important values
### constants -- THESE SHOULD BE IMPORTED FROM ASTROPY!


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
	#single_point_SNRs = np.sqrt((depths/hires_precision)**2)
	single_point_SNRs = depths / hires_precision ### equivalent to above, but not stupid.

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