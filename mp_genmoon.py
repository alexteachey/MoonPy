from __future__ import division
from pyluna import *
import matplotlib.pyplot as plt 


rhoSun = 1408.0
rhoJup = 1326.2
rhoNep = 1638.0
rhoEarth = 5515.0
RSun = 6.960e8
RJup = 69910558.8
RNep = 24621430.8
REarth = 6370400.0
q1Sun = 0.432984564534
q2Sun = 0.376524074763


class Moonpy_moon(object):
	

	def __init__(self, star_params=None, planet_params=None, sat_params=None, prompts='n', randomize='n', nmoons=1):
		self.nmoons = nmoons
		### either star_params, planet_params and sat_params are defined, OR the others must be.

		### STAR INPUTS
		if (star_params == None) and (prompts == 'n'):
			### use solar system properties
			self.rhostar, self.q1, self.q2 = rhoSun, q1Sun, q2Sun
			self.star_params = {'rhostar':self.rhostar, 'q1':self.q1, 'q2':self.q2}

		elif (star_params == None) and (prompts == 'y'):
			### prompt the user for inputs!
			self.rhostar = float(input('stellar density (solar is 1408 g/cm^3): '))
			self.q1 = float(input('q1 (solar is 0.43298): '))
			self.q2 = float(input('q2 (slar is 0.3765): '))
			self.star_params = {'rhostar':self.rhostar, 'q1':self.q1, 'q2':self.q2}

		else:
			self.star_params = star_params ### dictionary
			self.rhostar, self.q1, self.q2 = self.star_params['rhostar'], self.star_params['q1'], self.star_params['q2']


		### PLANET INPUTS 
		if (planet_params == None) and (prompts == 'n'):
			### make it Jupiter, with an earth year orbit
			self.rprstar = RJup / RSun
			self.bplan = 0
			self.Pplan = 365.25
			self.tau0 = 100 ### arbitrary
			self.rhoplan = rhoJup 
			self.planet_params = {'rprstar':self.rprstar, 'bplan':self.bplan, 'Pplan':self.Pplan, 'tau0':self.tau0, 'rhoplan':self.rhoplan}

		elif (planet_params == None) and (prompts == 'y'):
			self.rprstar = float(input('rprstar: '))
			self.bplan = float(input('impact parameter: '))
			self.Pplan = float(input('planet period [days]: '))
			self.tau0 = float(input('reference transit time: '))
			self.rhoplan = float(input('planet density (Earth = 5515 g/cm^3, Jupiter = 1326.2): '))
			self.planet_params = {'rprstar':self.rprstar, 'bplan':self.bplan, 'Pplan':self.Pplan, 'tau0':self.tau0, 'rhoplan':self.rhoplan}			

		else:
			self.planet_params = planet_params ### dictionary
			self.rprstar = self.planet_params['rprstar']
			self.bplan = self.planet_params['bplan']
			self.Pplan = self.planet_params['Pplan']
			self.tau0 = self.planet_params['tau0']
			self.rhoplan = self.planet_params['rhoplan']



		### SATELLITE INPUTS: 
		if (sat_params == None) and (prompts == 'n'):
			self.phi_sat = np.random.choice(np.linspace(0,2*np.pi,1000), size=self.nmoons)
			self.cosi_sat = np.random.choice(np.linspace(0.0,0.0,1000), size=self.nmoons) ### make them all edge on!
			self.omega_sat = np.random.choice(np.linspace(0,2*np.pi,1000), size=self.nmoons)

			if self.nmoons == 1:
			### make it Callisto!
				self.asp = 25.93 ### a_sp of Callisto
				self.msp = 5.667e-5 ### mass of Callisto / mass of Jupiter
				self.rsp = 0.034477
				self.phi_sat, self.cosi_sat, self.omega_sat = self.phi_sat[0], self.cosi_sat[0], self.omega_sat[0]
			else:
				self.asp = np.random.choice(np.linspace(5,100,1000), size=self.nmoons)
				self.msp = np.random.choice(np.linspace(0.0,0.0,10), size=self.nmoons)
				self.rsp = np.random.choice(np.linspace(0.01,0.3,100), size=self.nmoons)

			self.sat_params = {'asp':self.asp, 'phi_sat':self.phi_sat, 'cosi_sat':self.cosi_sat, 'omega_sat':self.omega_sat, 'msp':self.msp, 'rsp':self.rsp}

		elif (sat_params == None) and (prompts == 'y'):
			print('# of moons = ', self.nmoons)

			asps = []
			rsps = []
			cosi_sats = []
			phi_sats = []
			omega_sats = []

			for moon_idx in np.arange(1,self.nmoons+1,1):
				asp = float(input('SEMI-MAJOR AXIS a_sp for moon # '+str(moon_idx)+': '))
				asps.append(asp)
				rsp = float(input('RADIUS RATIO r_sp for moon # '+str(moon_idx)+': '))
				rsps.append(rsp)
				cosi_sat = input('COSI_SAT for moon # '+str(moon_idx)+' [-1 to 3, or press enter for 0]: ')
				if cosi_sat == '':
					cosi_sat = 0.0
				else:
					cosi_sat = float(cosi_sat)
				cosi_sats.append(cosi_sat)
				phi_sat = input("PHI_SAT for moon # "+str(moon_idx)+" (or press enter for random): ")
				if phi_sat == '':
					phi_sat = np.random.choice(np.linspace(0,2*np.pi,1000))
				else:
					phi_sat = float(phi_sat)
				phi_sats.append(phi_sat)
				omega_sat = input("OMEGA_SAT for moon# "+str(moon_idx)+" (or press enter for random): ")
				if omega_sat == '':
					omega_sat = np.random.choice(np.linspace(0,2*np.pi,1000))
				else:
					omega_sat = float(omega_sat)
				omega_sats.append(omega_sat)
				print(" ")

			self.asp = np.array(asps).astype(float)
			self.rsp = np.array(rsps).astype(float)
			self.cosi_sat = np.array(cosi_sats).astype(float)
			self.phi_sat = np.array(phi_sats).astype(float)
			self.omega_sat = np.array(omega_sats).astype(float)

			if self.nmoons == 1:
				self.asp = self.asp[0]
				self.rsp = self.rsp[0]
				self.cosi_sat = self.cosi_sat[0]
				self.phi_sat = self.phi_sat[0]
				self.omega_sat = self.omega_sat[0]
				self.msp = float(input("Mass ratio between secondary and primary (<1): "))
			
			else:
				### set the mass ratios to zero! so you don't have to deal with TTVs. CHANGE THIS.
				self.msp = np.linspace(0.0,0.0,self.nmoons)

			self.sat_params = {'asp':self.asp, 'phi_sat':self.phi_sat, 'cosi_sat':self.cosi_sat, 'omega_sat':self.omega_sat, 'msp':self.msp, 'rsp':self.rsp}

		else:
			self.sat_params = sat_params ### dictionary
			self.asp = self.sat_params['asp']
			self.phi_sat = self.sat_params['phi_sat']
			self.cosi_sat = self.sat_params['cosi_sat']
			self.omega_sat = self.sat_params['omega_sat']
			self.msp = self.sat_params['msp']
			self.rsp = self.sat_params['rsp']



	def gen_transit(self, tau=None, window=None, cadence_minutes=None, prompts='n', ntransits=1, show_plots='n', ppm=0.0, model='M'):
		if (tau == None) and (prompts == 'n'):
			tau = self.tau0
		elif (tau == None) and (prompts == 'y'):
			tau = float(input('What is the reference time for the transit? (tau0 = '+str(self.tau0)+'): '))

		if (window == None) and (prompts == 'n'):
			window = 5
		elif (window == None) and (prompts == 'y'):
			window = float(input('What is the time window on either side of the transit? [days]: '))

		if (cadence_minutes == None) and (prompts == 'n'):
			cadence_minutes = 29.42 ### native Kepler cadence
		elif (cadence_minutes == None) and (prompts == 'y'):
			cadence_minutes = float(input('What is the cadence? [minutes]: '))


		### generate a single transit centered around some tau
		cadence_days = cadence_minutes / (60 * 24)
		if ntransits == 1:
			all_times = np.arange(tau - window, tau + window + cadence_days, cadence_days)
		else:
			all_times = []
			for transitnum in np.arange(0,ntransits,1):
				tmid = tau+(self.Pplan*transitnum)
				all_times.append(np.arange(tmid-window, tmid+window+cadence_days, cadence_days))
			all_times = np.hstack(all_times)

		if model == "M":
			nparamorig = 14
			nparam = 14
			nvars = 14 ### fitting all the parameters!
		elif model == "P":
			#nparamorig = 8  
			#nparam = 8
			nparamorig = 14 ### all these inputs must still be present, even if some of them are fixed at zero!
			nparam = 14
			nvars = 8  ### RpRstar, rhostar, bplan, Pplan, tau0, q1, q2, rhoplan
		elif model == 'T':
			#nparamorig = 8 ### RpRstar, rhostar, bplan, Pplan, tau0, q1, q2, rhoplan 
			nparamorig = 14
			nparam = nparamorig + (ntaus-1) ### tau0 is a STANDARD nparamorig input... every additional tau needs to be counted.
			nvars = 8 + (ntaus-1) ### standard P model variables plus all additional taus.
		elif model == 'Z':
			nparam = 14
			nparamorig = 14
			nvars = 13 ### not fitting Rsat/Rp 

		### calculate ntaus
		ntaus = 0
		first_tau = self.tau0
		while first_tau <= np.nanmin(all_times):
			first_tau = first_tau + self.Pplan
		last_tau = self.tau0
		while last_tau <= np.nanmax(all_times):
			last_tau = last_tau + self.Pplan
		if last_tau > np.nanmax(all_times):
			last_tau = last_tau - self.Pplan
		
		all_taus = np.arange(first_tau, last_tau+self.Pplan, self.Pplan)
		ntaus = len(all_taus)


		prepare_files(all_times, ntaus, nparam, nparamorig)

		if self.nmoons == 1:
			output_times, output_fluxes = run_LUNA(all_times, self.rprstar, self.rhostar, self.bplan, self.Pplan, self.tau0, self.q1, self.q2, self.rhoplan, self.asp, self.phi_sat, self.cosi_sat, self.omega_sat, self.msp, self.rsp, model="M", cadence_minutes=cadence_minutes, print_params='y')
		
		else:
			### generate a planet only light curve first!
			planet_only_times, planet_only_fluxes = run_LUNA(all_times, self.rprstar, self.rhostar, self.bplan, self.Pplan, self.tau0, self.q1, self.q2, self.rhoplan, 1000, 0, 0, 0, 0, 0, model='M', cadence_minutes=cadence_minutes, print_params='y')
			planet_only_missing_fluxes = 1 - planet_only_fluxes

			if show_plots=='y':
				fig, (ax1, ax2) = plt.subplots(2, sharex=True)
				ax1.scatter(planet_only_times, planet_only_fluxes, facecolor='LightCoral', edgecolor='k', s=10)
				ax1.set_ylabel('Flux')
				ax2.scatter(planet_only_times, planet_only_missing_fluxes, facecolor='LightCoral', edgecolor='k', s=10)
				ax2.set_ylabel("Missing Flux")
				plt.show()

			for moon_idx in np.arange(0,self.nmoons,1):
				moon_and_planet_times, moon_and_planet_fluxes = run_LUNA(all_times, self.rprstar, self.rhostar, self.bplan, self.Pplan, self.tau0, self.q1, self.q2, self.rhoplan, self.asp[moon_idx], self.phi_sat[moon_idx], self.cosi_sat[moon_idx], self.omega_sat[moon_idx], self.msp[moon_idx], self.rsp[moon_idx], model="M", cadence_minutes=cadence_minutes, print_params='y')
				moon_only_fluxes = moon_and_planet_fluxes / planet_only_fluxes
				moon_only_missing_fluxes = 1 - moon_only_fluxes 

				if show_plots=='y':
					fig,(ax1, ax2, ax3) = plt.subplots(3, sharex=True)
					ax1.scatter(moon_and_planet_times, moon_and_planet_fluxes, facecolor='LightCoral', edgecolor='k', s=10)
					ax1.set_ylabel('moon+planet')
					ax1.set_ylim(np.nanmin(moon_and_planet_fluxes), np.nanmax(moon_and_planet_fluxes))
					ax2.scatter(moon_and_planet_times, moon_only_fluxes, facecolor='LightCoral', edgecolor='k', s=10)
					ax2.set_ylabel('moon only')
					ax2.set_ylim(np.nanmin(moon_only_fluxes), np.nanmax(moon_only_fluxes))
					ax3.scatter(moon_and_planet_times, moon_only_missing_fluxes, facecolor='LightCoral', edgecolor='k', s=10)
					ax3.set_ylabel('missing moon flux')
					ax3.set_ylim(np.nanmin(moon_only_missing_fluxes), np.nanmax(moon_only_missing_fluxes))
					plt.show()
			
				if moon_idx == 0:
					moon_only_flux_stack = moon_only_fluxes
					moon_only_missing_flux_stack = moon_only_missing_fluxes
				else:
					moon_only_flux_stack = np.vstack((moon_only_flux_stack, moon_only_fluxes))
					moon_only_missing_flux_stack = np.vstack((moon_only_missing_flux_stack, moon_only_missing_fluxes))

			missing_fluxes = planet_only_missing_fluxes + np.nansum(moon_only_missing_flux_stack, axis=0)
			final_fluxes = 1 - missing_fluxes

			if ppm != 0.0:
				final_fluxes = np.random.normal(loc=final_fluxes, scale=ppm*1e-6)


			if show_plots=='y':
				plt.scatter(planet_only_times, final_fluxes, facecolor='LightCoral', edgecolor='k')
				plt.xlabel('Time')
				plt.ylabel('Flux')
				plt.show()

			output_times = planet_only_times 
			output_fluxes = final_fluxes 
			### verify this is working as intended


		self.times, self.fluxes = output_times, output_fluxes

	def plot_moon(self):
		try:
			print(self.times)
		except:
			self.gen_transit()

		plt.scatter(self.times, self.fluxes, facecolor='LightCoral', edgecolor='k', s=20)
		plt.xlabel('Time')
		plt.ylabel('Flux')
		plt.show()





		