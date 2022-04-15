from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pandas
import os
import traceback
from astropy.io import fits
try:
	import exoplanet as xo
except:
	print('COULD NOT IMPORT the exoplanet package.')
	print("type 'conda install -c conda-forge exoplanet' to install.")
	print(' ')
import astropy.units as u 
try:
	import pymc3 as pm
except:
	print("COULD NOT IMPORT pymc3.")
	print("type 'conda install -c conda-forge pymc3' to install.")
	print(' ')
try:
	import pymc3_ext as pmx
except:
	print('COULD NOT IMPORT pymc3_ext')
	print("type 'conda install -c conda-forge pymc3_ext' to install.")
	print(' ')

try:
	import arviz as az
except:
	print('COULD NOT IMPORT arviz.')
	print("type 'conda install -c conda-forge arviz=0.11.0' to install.")
	print(' ')
try:
	import corner
except:
	print('COULD NOT IMPORT corner.')
	print("type conda install -c astropy corner' to install.")
	print(' ')

try:
	from celerite2.theano import terms, GaussianProcess
except:
	print('COULD NOT IMPORT celerite2 modules.')
	print("type 'conda install -c conda-forge celerite2' to install.")
	print(' ')

try:
	import aesara_theano_fallback.tensor as tt
except:
	print("COULD NOT IMPORT aesara_theano_fallback.tensor")
	print("type 'conda install -c conda-forge aesara-theano-fallback' to install.")
	print(' ')

import platform 
from matplotlib import rcParams
from moonpy import * 
import pickle

#rcParams['font.family'] = 'serif'


#### THIS CODE WILL TAKE CHETAN'S VETTING, GRAB THE LIGHT CURVES, PREPARE THEM FOR THE EXOPLANET RUN.

#show_test_plots = input('Show test plots? y/n: ')
show_test_plots = 'n'



ncores = 2
nchains = 2
ndraws = 1000
nsamples = ndraws * nchains

#skip_completed = input("Do you want to skip planets you've already run? y/n: ")

moonpydir = os.path.realpath(__file__)
moonpydir = moonpydir[:moonpydir.find('/_mp_planet_fitter.py')]

hostname = socket.gethostname()
if ('tethys' in hostname) and ('sinica' in hostname):
	#moonpydir = '/data/tethys/Documents/Software/MoonPy'
	central_data_dir = '/data/tethys/Documents/Central_Data'
elif ('Alexs-MacBook') in hostname:
	#moonpydir = '/Users/hal9000/Documents/Software/MoonPy'
	central_data_dir = '/Users/hal9000/Documents/Central_Data'
elif 'umbriel' in hostname:
	#moonpydir = '/home/cal/ateachey/Documents/MoonPy'
	central_data_dir = '/home/cal/ateachey/Documents/Central_Data'
else:
	#moonpydir = input('Please specify the MoonPy directory (or hard-code this into moonpy.py): ')
	#central_data_dir = input("Please specify a 'central data' directory (or hard-code this into moonpy.py): ")
	### store central_data within MoonPy directory
	if os.path.exists(moonpydir+'/Central_Data'):
		pass
	else:
		os.system('mkdir '+moonpydir+'/Central_Data')
	central_data_dir = moonpydir+'/Central_Data'



if os.path.exists(moonpydir+'/fitting_errors.txt') == False:
	#### make it!
	os.system('touch '+moonpydir+'/fitting_errors.txt')
else:
	pass



#kicfile = pandas.read_csv('/home/amteachey/Documents/Projects/Unistellar_transit_obs/two_transits.csv', encoding='unicode_escape')
#kicdict = {}

#for col in kicfile.columns:
#	kicdict[col] = np.array(kicfile[col])

#kics = []
#for kic in kicdict['KIC']:
#	kics.append('KIC'+str(kic)[:-3])
#kics = np.array(kics)


#for nvpi,vpi in enumerate(vetting_priority_idxs):
#for nkic,kic in enumerate(kics):


def run_planet_fit(self, period=None, tau0=None, tdur_hours=None, smass=None, smass_err=None, show_plots=True, use_mp_detrend=False, restrict_data=True, fit_neighbors=False, savepath=None):

	print('fit_neighbors: ', fit_neighbors)

	number_of_neighbors = len(self.neighbor_dict.keys())
	if fit_neighbors == True:
		total_number_of_planets = number_of_neighbors+1 
	elif fit_neighbors == False:
		total_number_of_planets = 1 


	if show_plots == True:
		keep_showing = 'y'
	else:
		keep_showing = 'n'

	if type(savepath) == type(None):
		savepath = self.savepath 


	if type(period) != type(None):
		#### update it
		self.period = period 

	if type(tau0) != type(None):
		self.tau0 = tau0 

	if type(tdur_hours) != type(None):
		### update it
		self.duration_hours = tdur_hours
		self.duration_days = self.duration_hours / 24 

	if type(smass) != type(None):
		#### update it
		self.smass = smass


	if type(smass_err) != type(None):
		#### update it
		self.smass_err = smass_err
	else:
		self.smass_err = 0.05*self.smass 


	try:
		print('self.period = ', self.period)
	except:
		manual_period_entry = input('enter period: ')
		self.period = float(manual_period_entry)

	try:
		print('self.tau0 = ', self.tau0)
	except:
		manual_tau0_entry = input('enter tau0: ')
		self.tau0 = float(manual_tau0_entry)

	try:
		print('self.duration_hours = ', self.duration_hours)
	except:
		manual_duration_hours_entry = input('Enter transit duration, in hours: ')
		self.duration_hours = float(manual_duration_hours_entry)

	try:
		print('self.smass = ', self.smass)
		if np.isfinite(self.smass) == False:
			manual_smass_entry = input('Enter stellar mass, in solar units: ')
			self.smass = float(manual_smass_entry)

	except:
		manual_smass_entry = input('Enter stellar mass, in solar units: ')
		self.smass = float(manual_smass_entry)
		

	try:
		print('self.smass_err = ', self.smass_err)
		if np.isfinite(self.smass_err) == False:
			manual_smass_err_entry = input('Enter stellar mass err, in solar units: ')
			self.smass_err = float(manual_smass_entry)
	except:
		manual_smass_err_entry = input('Enter stellar mass err, in solar units: ')
		self.smass_err = float(manual_smass_entry)
		




	try:

		#duration = tdur_days  #### I think we should use this...
		#tau0 = tau0 


		if use_mp_detrend == True:
			times, fluxes, errors, fluxes_detrend, errors_detrend, flags = self.times, self.fluxes, self.errors, self.fluxes_detrend, self.errors_detrend, self.flags

		elif use_mp_detrend == False:
			times, fluxes, errors, flags = self.times, self.fluxes, self.errors, self.flags
			fluxes_detrend, errors_detrend = [], []

			for f,e in zip(fluxes, errors):
				#### have to do this on a quarter-by-quarter basis 
				fd = f/np.nanmedian(f)
				fluxes_detrend.append(fd)

				ed = e/f 
				errors_detrend.append(ed)

			fluxes_detrend, errors_detrend = np.array(fluxes_detrend), np.array(errors_detrend)
			self.fluxes_detrend, self.errors_detrend = fluxes_detrend, errors_detrend 


		taus = [self.tau0]
		while taus[-1] < np.nanmax(np.concatenate(times)):
			taus.append(taus[-1] + self.period)
		taus = np.array(taus)
		self.taus = taus 	



		cctimes, ccsap, ccsap_err, cc_fluxes_detrend, cc_errors_detrend, ccflags = np.concatenate(times), np.concatenate(fluxes), np.concatenate(errors), np.concatenate(fluxes_detrend), np.concatenate(errors_detrend), np.concatenate(flags)
		good_flags = np.where(ccflags == 0)[0]


		#### find nearby fluxes -- don't want to fit everything!!!!! just use +/- 3 days on either side of a transit
		if restrict_data == True:
			if fit_neighbors == False:
				near_transit_idxs = []
				for ncctime, cctime in enumerate(cctimes):
					if np.any(np.abs(cctime - taus) < 10):
						near_transit_idxs.append(ncctime)
				near_transit_idxs = np.array(near_transit_idxs)

			elif fit_neighbors == True:
				#### just need to the whole light curve
				near_transit_idxs = np.arange(0,len(cctimes),1)


			print('len(cctimes) = ', len(cctimes))
			print('len(near_transit_idxs) = ', len(near_transit_idxs))
			continue_query = input('Do you wish to continue? y/n: ')
			if continue_query != 'y':
				raise Exception('you opted not to continue.')

			##### update !
			cctimes, ccsap, ccsap_err, cc_fluxes_detrend, cc_errors_detrend, ccflags = cctimes[near_transit_idxs], ccsap[near_transit_idxs], ccsap_err[near_transit_idxs], cc_fluxes_detrend[near_transit_idxs], cc_errors_detrend[near_transit_idxs], ccflags[near_transit_idxs] 




		#### now compute some expected values from the BLS results
		first_transit = self.tau0
		if first_transit > np.nanmin(cctimes) + self.period:
			while first_transit > np.nanmin(cctimes) + self.period:
				first_transit = first_transit - self.period 
		
		elif first_transit < np.nanmin(cctimes):
			while first_transit < np.nanmin(cctimes):
				first_transit = first_transit + self.period 

		else:
			pass

		transit_times = np.arange(first_transit, np.nanmax(cctimes), self.period)
		print("transit_times = ", transit_times)


		##### MASK OUT THE TRANSITS FOR GP FITTING.
		transit_idxs = []
		out_of_transit_idxs = []
		for ncct,cct in enumerate(cctimes):
			if np.any(np.abs(cct - transit_times) < 0.5*self.duration_days):
				#### means it's in transit
				transit_idxs.append(ncct)
			else:
				out_of_transit_idxs.append(ncct)

		transit_idxs, out_of_transit_idxs = np.array(transit_idxs), np.array(out_of_transit_idxs)




		#### CREATE A PHASE-FOLD
		fold_times = []
		for ctt in cctimes:
			foldt = ctt
			if foldt > first_transit + (0.5*self.period):
				#### subtract off a period until you get there
				while foldt > first_transit + (0.5*self.period):
					foldt = foldt - self.period 
			elif foldt < first_transit - (0.5*self.period):
				#### add a period until you get there
				while foldt < first_transit - (0.5*self.period):
					foldt = foldt+ self.period 
			fold_times.append(foldt)
		fold_times = np.array(fold_times, dtype=object)
		fold_times = fold_times - first_transit	



		#### now plot them all up and see how well BLS is doing
		if show_test_plots == 'y':

			fig, ax = plt.subplots(2)

			ax[0].scatter(cctimes, cc_fluxes_detrend, facecolor='LightCoral', edgecolor='k', s=20, alpha=0.7)

			for btt in transit_times:
				#### mark with a vertical line
				ax[0].plot(np.linspace(btt,btt,100), np.linspace(0.95*np.nanmin(cc_fluxes_detrend), 1.05*np.nanmax(cc_fluxes_detrend), 100), color='red', alpha=0.7, linewidth=1)
				#plt.plot(np.linspace(btt-0.5*duration,btt-0.5*duration,100), np.linspace(0.95*np.nanmin(cc_fluxes_detrend), 1.05*np.nanmax(cc_fluxes_detrend), 100), color='red', alpha=0.8, linestyle='--')
				#plt.plot(np.linspace(btt+0.5*duration,btt+0.5*duration,100), np.linspace(0.95*np.nanmin(cc_fluxes_detrend), 1.05*np.nanmax(cc_fluxes_detrend), 100), color='red', alpha=0.8, linestyle='--')			

			ax[1].scatter(fold_times, ccksp, facecolor='LightCoral', edgecolor='k', s=20, alpha=0.7)

			ax[0].set_title(kic)
			#ax[1].set_xlabel('BTJD')
			ax[0].set_ylabel('FLUX')
			ax[1].set_ylabel('FLUX')
			plt.show()




		#### NOW THE TRANSIT FITTING!
		#### The transit model in PyMC3

		#### MY TRY
		if fit_neighbors == False:
			planet_names = self.target
			periods = self.period #### my guess
			t0s = first_transit
			rprstars = self.rprstar


		elif fit_neighbors == True:
			planet_names = [self.target]
			periods = [self.period]
			t0s = [self.tau0]
			rprstars = [self.rprstar]

			for key in self.neighbor_dict.keys():
				planet_names.append(self.neighbor_dict[key].target)
				periods.append(self.neighbor_dict[key].period)
				t0s.append(self.neighbor_dict[key].tau0)
				rprstars.append(self.neighbor_dict[key].rprstar)

			planet_names = np.array(planet_names)
			periods = np.array(periods)
			t0s = np.array(t0s)
			rprstars = np.array(rprstars)

			assert len(periods) == total_number_of_planets
			assert len(t0s) == total_number_of_planets 
			assert len(planet_names) == total_number_of_planets
			assert len(rprstars) == total_number_of_planets 



		model_times = cctimes
		yvals = cc_fluxes_detrend
		yerr = cc_errors_detrend
		

		"""
		YOU SHOULD REALLY FOLLOW THIS PAGE -- TO INCORPORATE THE GPs!!!
		https://gallery.exoplanet.codes/tutorials/tess/#the-transit-model-in-pymc3

		"""


		with pm.Model() as model: #### this is a PYMC3 MODEL! 

			phase_lc = np.linspace(-self.period/2, self.period/2, 10000) #### what is this?!


			#### COLLECTING ALL THE PRIORS HERE.

			#### LIGHT CURVE PRIORS 
			mean = pm.Normal("mean", mu=1.0, sd=np.nanmedian(yerr))	
			ldcs = xo.distributions.QuadLimbDark("ldcs", testval=(0.3, 0.2)) #### making it a tuple -- I guess this is what they want?
			
			#### STELLAR PRIORS
			star = xo.LimbDarkLightCurve(ldcs[0], ldcs[1])
			BoundedNormal = pm.Bound(pm.Normal, lower=0.0)
			try:
				m_star = BoundedNormal("m_star", mu=self.smass, sd=np.nanmean(np.abs(self.smass_err)))
			except:
				m_star = BoundedNormal("m_star", mu=self.smass, sd=0.05*self.smass)
			r_star = BoundedNormal("r_star", mu=self.rstar_rsol, sd=0.05*self.rstar_rsol)


			#### PLANET PRIORS 
			t0 = pm.Normal("t0", mu=t0s, shape=total_number_of_planets, sd=1.0)
			log10Period = pm.Normal("log10Period", mu=np.log10(periods), shape=total_number_of_planets, sd=0.01)
			period = pm.Deterministic("period", 10**log10Period)	
			impact = pm.Uniform("impact", lower=0, upper=1, shape=total_number_of_planets) 
			log_depth = pm.Normal("log_depth", mu=np.log(np.nanmin(yvals)), shape=total_number_of_planets, sigma=2.0) ### why sigma=2?
			try:
				rprstar = pm.Uniform('rprstar', lower=0, upper=1, shape=total_number_of_planets, testval=rprstars)
			except:
				rprstar = pm.Uniform('rprstar', lower=0, upper=1, shape=total_number_of_planets, testval=0.05)
			rplan = pm.Deterministic("rplan", rprstars * r_star) 

			#### THIS MAY NEED SOME ADJUSTING!!!
			if fit_neighbors == True:
				ecs = pmx.UnitDisk('ecs', testval=(0.01 * np.ones((2,total_number_of_planets))), shape=(2,total_number_of_planets)) #### making it a tuple... I think this is what they want.			
				ecc = pm.Deterministic("ecc", tt.sum(ecs ** 2, axis=0))		
			
			elif fit_neighbors == False:
				ecs = pmx.UnitDisk('ecs', testval=(0.01, 0.0)) #### ORIGINAL -- WORKS
				ecc = pm.Deterministic("ecc", tt.sum(ecs ** 2))		


			omega = pm.Deterministic("omega", tt.arctan2(ecs[1], ecs[0])) #### what is this??!
			xo.eccentricity.kipping13("ecc_prior", fixed=True, observed=ecc, shape=total_number_of_planets)
			

			#Transit jitter & GP parameters
			log_sigma_lc = pm.Normal("log_sigma_lc", mu=np.log(np.nanstd(yvals)), sd=10) #### why sd=10?
			log_rho_gp = pm.Normal("log_rho_gp", mu=0, sd=10) #### what is this mean and sd?!?!
			log_sigma_gp = pm.Normal("log_sigma_gp", mu=np.log(np.nanstd(yvals)), sd=10)  #### 

			# Set up a Keplerian orbit for the planets -- models the system
			orbit = xo.orbits.KeplerianOrbit(period=period, t0=t0, b=impact, ecc=ecc, omega=omega, m_star=m_star, r_star=r_star)

			#### creates an observation from the model
			light_curves = star.get_light_curve(orbit=orbit, r=rplan, t=model_times)	
			#### trying a MASK!
			#light_curves = star.get_light_curve(orbit=orbit, r=rplan, t=model_times)			
			light_curve = tt.sum(light_curves, axis=-1) + mean 

			residuals = yvals - light_curve 
			#### TRYING A MASK!
			#residuals = yvals - light_curve 


			### GP MODEL FOR THE LIGHT CURVE
			kernel = terms.SHOTerm(sigma=tt.exp(log_sigma_gp), rho=tt.exp(log_rho_gp), Q=1 / np.sqrt(2),)
			#gp = GaussianProcess(kernel, t=model_times, yerr=tt.exp(log_sigma_lc))
			#### trying a MASK! 
			gp = GaussianProcess(kernel, t=model_times, yerr=tt.exp(log_sigma_lc))
			gp.marginal("gp", observed=residuals)
			

			# Here we track the value of the model light curve for plotting
			# purposes
			pm.Deterministic("light_curves", light_curves)
			if fit_neighbors == False:
				pm.Deterministic("lc_pred", 1e6 * star.get_light_curve(orbit=orbit, r=rplan, t=t0 + phase_lc)[...,0]) #### WHAT IS ALL THIS?!
			elif fit_neighbors == True:
				### see this page for the syntax example: https://gallery.exoplanet.codes/tutorials/joint/ 
				pm.Deterministic(
					"lc_pred", 
					1e6 * tt.stack(
						[
						star.get_light_curve(
							orbit=orbit, r=rplan, t=t0[n] + phase_lc
							)[...,n]
							for n in range(total_number_of_planets)
							],
							axis=-1,
					),
				) 




			# ******************************************************************* #
			# On the folowing lines, we simulate the dataset that we will fit     #
			#                                                                     #
			# NOTE: if you are fitting real data, you shouldn't include this line #
			#       because you already have data!                                #
			# ******************************************************************* #
			#y = pmx.eval_in_model(light_curve)
			#y += yerr * np.random.randn(len(y))
			# ******************************************************************* #
			# End of fake data creation; you want to include the following lines  #
			# ******************************************************************* #



			# The likelihood function assuming known Gaussian uncertainty
			#pm.Normal("obs", mu=light_curve, sd=yerr, observed=y) ### UNNECESSARY?!

			# Fit for the maximum a posteriori parameters given the simuated
			# dataset --- HOW DOES THIS WORK
			#map_soln = pmx.optimize(start=model.test_point)

			map_soln = pmx.optimize(start=model.test_point, vars=[log_sigma_lc, log_sigma_gp, log_rho_gp])
			map_soln = pmx.optimize(start=map_soln, vars=[log_depth])
			map_soln = pmx.optimize(start=map_soln, vars=[impact])
			map_soln = pmx.optimize(start=map_soln, vars=[log10Period, t0])
			map_soln = pmx.optimize(start=map_soln, vars=(log10Period, t0))
			map_soln = pmx.optimize(start=map_soln, vars=[ldcs])
			map_soln = pmx.optimize(start=map_soln, vars=[log_depth]) ### huh? why again?
			map_soln = pmx.optimize(start=map_soln, vars=[impact]) #### why again?
			map_soln = pmx.optimize(start=map_soln, vars=[ecs])
			map_soln = pmx.optimize(start=map_soln, vars=[mean])
			map_soln = pmx.optimize(start=map_soln, vars=[log_sigma_lc, log_sigma_gp, log_rho_gp]) #huh?
			map_soln = pmx.optimize(start=map_soln) 
			

			extras = dict(zip(["light_curves", "gp_pred"], pmx.eval_in_model([light_curves, gp.predict(residuals)], map_soln),))


		self.exoplanet_model = model 
		self.exoplanet_map_soln = map_soln 
		self.exoplanet_extras = extras 


		model0, map_soln0, extras0 = model, map_soln, extras #### FOR DEBUGGING!

		#### PLOT THE MAXIMUM A POSTERIORI SOLUTION FROM THE FITTING ABOVE
		fig, axes = plt.subplots(3, 1, figsize=(10,7), sharex=True)

		ax = axes[0]
		ax.scatter(model_times, yvals*1e6, facecolor='LightCoral', edgecolor='k', marker='o', s=20, alpha=0.5, label="data", zorder=0)	
		gp_mod = extras["gp_pred"] + map_soln["mean"]
		ax.plot(model_times, gp_mod*1e6, color="C2", label="gp model")
		ax.legend(fontsize=10, loc='upper right')
		ax.set_ylabel("relative flux [ppm]")	

		ax = axes[1]
		ax.scatter(model_times, (yvals - gp_mod)*1e6, facecolor='LightCoral', edgecolor='k', marker='o', s=20, alpha=0.5, zorder=0)
		for i in np.arange(0,total_number_of_planets,1):
			mod = extras["light_curves"][:,i]
			scatter_ppm = np.nanstd((yvals - mod)*1e6)
			ax.plot(model_times, mod*1e6, label=planet_names[i], zorder=1)		
		
		ax.set_ylim(np.nanmin((yvals-gp_mod)*1e6) - 5*scatter_ppm, 5*scatter_ppm)
		ax.legend(fontsize=10, loc='upper right')
		ax.set_ylabel("de-trended flux [ppm]")


		ax = axes[2]
		mod = gp_mod + np.sum(extras["light_curves"], axis=-1)
		ax.scatter(model_times, (yvals - mod)*1e6, facecolor='LightCoral', edgecolor='k', marker='o', s=20, alpha=0.5, zorder=0)
		ax.axhline(0, color="#aaaaaa", lw=1)
		ax.set_ylabel("residuals [ppm]")
		ax.set_xlim(model_times.min(), model_times.max())
		ax.set_ylim(-5*scatter_ppm, 5*scatter_ppm)
		ax.set_xlabel("time [days]")



		#plt.savefig(kicpath+'/'+kic+'_gp_detrend_and_residuals.pdf', dpi=300)
		plt.savefig(savepath+'/'+self.target+'_gp_detrend_and_residuals.pdf', dpi=300)
		if keep_showing == 'y':
			plt.show()
		plt.close()



		with model:
			trace = pm.sample(
				tune=1500,
				draws=ndraws,
				start=map_soln,
				cores=ncores,
				chains=nchains,
				target_accept=0.95,
				return_inferencedata=True,
				init='adapt_full',
			)



		az.summary(
			trace,
			var_names=[
				"omega",
				"ecc",
				"rplan",
				"impact",
				"t0",
				"period",
				"r_star",
				"m_star",
				"ldcs",
				"mean",
			],
		)



		#### PLOT FROM DRAWS
		flat_samps = trace.posterior.stack(sample=("chain", "draw"))

		with open(savepath+'/'+self.target+'_flat_samps.pkl', 'wb') as handle:
			pickle.dump(flat_samps, handle)

		self.flat_samps = flat_samps #### so it can be accessed later

		# Compute the GP prediction
		gp_mod = extras["gp_pred"] + map_soln["mean"]  # np.median(
		#     flat_samps["gp_pred"].values + flat_samps["mean"].values[None, :], axis=-1
		# )

		# Get the posterior median orbital parameters
		post_period = np.median(flat_samps["period"])
		post_t0 = np.median(flat_samps["t0"])

		# Plot the folded data
		x_fold = (model_times - post_t0 + (0.5 * post_period)) % post_period - (0.5 * post_period)
		yvals_ppm = (yvals - gp_mod) * 1e6
		plt.scatter(x_fold, yvals_ppm, label="data", facecolor='LightCoral', edgecolor='k', alpha=0.2, s=10, zorder=-1000)

		# Plot the folded model
		pred = np.percentile(flat_samps["lc_pred"], [16, 50, 84], axis=-1) ##### THIS IS WHERE IT ALL GOES TO SHIT -- WHAT IS THE PROBLEM?!?!
		plt.plot(phase_lc, pred[1], color="C1", label="model")
		art = plt.fill_between(
			phase_lc, pred[0], pred[2], color="C1", alpha=0.5, zorder=1000
		)
		art.set_edgecolor("none")

		# Annotate the plot with the planet's period
		txt = "period = {0:.5f} +/- {1:.5f} d".format(
			np.mean(flat_samps["period"].values), np.std(flat_samps["period"].values)
		)
		plt.annotate(
			txt,
			(0, 0),
			xycoords="axes fraction",
			xytext=(5, 5),
			textcoords="offset points",
			ha="left",
			va="bottom",
			fontsize=12,
		)

		plt.legend(fontsize=10, loc=4)
		#plt.xlim(-0.5 * post_period, 0.5 * post_period)
		plt.xlabel("time since transit [days]")
		plt.ylabel("de-trended flux [ppm]")
		try:
			plt.xlim(-3*self.duration_days, 3*self.duration_days)
		except:
			pass
		#_ = plt.xlim(-0.15, 0.15)
		#_ = plt.ylim(np.nanmin(yvals_ppm)-scatter_ppm, scatter_ppm)

		plt.subplots_adjust(left=0.15, right=0.85, top=0.9, bottom=0.1)

		#plt.savefig(kicpath+'/'+kic+'_phasefold.pdf', dpi=300)
		plt.savefig(savepath+'/'+self.target+'_phasefold.pdf', dpi=300)
		if keep_showing == 'y':
			plt.show()
		plt.close()




		trace.posterior["r_earth"] = (
			trace.posterior["rplan"].coords,
			(trace.posterior["rplan"].values * u.R_sun).to(u.R_earth).value,
		)

		trace.posterior['ldc1'] = trace.posterior['ldcs'][:,:,0]
		trace.posterior['ldc2'] = trace.posterior['ldcs'][:,:,1]



		#### PHYSICAL CORNER PLOT
		try:
			cornerfig = corner.corner(
				trace,
				var_names=["period", "m_star", "r_star", "rprstar", "ldc1", "ldc2", "impact", "ecc", "omega"],
				labels=[
					r"$P$",
					r"$M_{*} \, [\odot]$",
					r"$R_{*} \, [\odot]$",
					r"$R_{P} / R_{*}$",
					r'$q_1$',
					r'$q_2$',
					r"$b$",
					r"$e$",
					r'$\omega$',
				],
			)


			#plt.savefig(kicpath+'/'+kic+'_physical_cornerplot.pdf', dpi=300)
			plt.savefig(savepath+'/'+self.target+'_physical_cornerplot.pdf', dpi=300)
			if keep_showing == 'y':
				plt.show()
			plt.close()

		except:
			print('COULD NOT PRODUCE THE CORNER PLOT. MAY NEED TO INSTALL CORNER.')



		##### MODEL CORNER PLOT (GP AND JITTER)
		#Transit jitter & GP parameters
		cornerfig = corner.corner(
			trace,
			var_names=["log_sigma_lc", "log_rho_gp", "log_sigma_gp"],
			labels=[
				r'$\sigma_{\mathrm{LC}}$',
				r'$\rho_{\mathrm{GP}}$',
				r'$\sigma_{\mathrm{GP}}$',
			],
		)

		#plt.savefig(kicpath+'/'+kic+'_jitter_cornerplot.pdf', dpi=300)
		plt.savefig(savepath+'/'+self.target+'_jitter_cornerplot.pdf', dpi=300)
		if keep_showing == 'y':
			plt.show()
		plt.close()



		#### all the parameters you've fit
		flat_samps_keys = ['mean','ldc1','ldc2','m_star','r_star','log_sigma_lc','log_rho_gp','log_sigma_gp','t0','log10Period','period','impact','log_depth','rprstar','rplan','ecs1','ecs2','ecc','omega']
		#posterior_file = open(kicpath+'/'+kic+'_posteriors.csv', mode='w')
		
		headers = []
		star_headers = ['mean','ldc1','ldc2','m_star','r_star','log_sigma_lc','log_rho_gp','log_sigma_gp'] #### stellar parameters, not planet-specific 



		for fsk in flat_samps_keys:
			if fsk in star_headers:
				headers += [fsk]

			else:
				### means it's a planet value -- do these individually!
				for i in np.arange(0,total_number_of_planets,1):
					#### append the planet number to the header name -- all in a row!
					headers += [fsk+'_'+str(i)]




		#### now we have a whole bunch of headers -- for the planets, they'll be aranged t0_0, t0_1, log10Period_0, log10Period_1, etc...

		posterior_file = open(savepath+'/'+self.target+'_posteriors.csv', mode='w')
		
		### write the header
		#for nfsk,fsk in enumerate(flat_samps_keys):
		for nheader,header in enumerate(headers):
			if header != headers[-1]:
				posterior_file.write(header+',')
			elif header == headers[-1]:
				posterior_file.write(header+'\n')


		continue_printing = 'y'


		for sample in np.arange(0,nsamples,1):
			##### for each draw of the posterior
			for nfsk,fsk in enumerate(flat_samps_keys):
				#try:
				#### pull out the value for that particular parameter

				if fsk in star_headers:
					#### means there is only a single value for the system, not planet-specific.
					if fsk == 'ldc1':
						sampval = np.array(flat_samps['ldcs'])[0][sample]
					elif fsk == 'ldc2':
						sampval = np.array(flat_samps['ldcs'])[1][sample]	

					else:
						sampval = np.array(flat_samps[fsk][sample])

					posterior_file.write(str(sampval)+',') #### all the star values come first!


				else:
					#### planet values
					if total_number_of_planets == 1:
						#### don't need the extra indexing
						if fsk == 'ecs1':
							sampval = np.array(flat_samps['ecs'])[0][sample]						
						elif fsk == 'ecs2':
							sampval = np.array(flat_samps['ecs'])[1][sample]

						else:
							try:
								#sampval = np.array(flat_samps[fsk])[sample][sample] #### there's only one!
								sampval = np.array(flat_samps[fsk])[0][sample]
							except:
								try:
									sampval = np.array(flat_samps[fsk])[sample]
								except:
									traceback.print_exc()
									print(' ')
									print(' ')
									print('fsk = ', fsk)
									print('sample = ', sample)
									print('len(flat_samps[fsk] = ', len(flat_samps[fsk]))
									posterior_file.close()
									print('Posteriors written to file.')

							#try:
							#	sampval = np.array(flat_samps[fsk])[0][sample]
							#except:
							#	sampval = np.array(flat_samps[fsk])[sample]


					else:
						for i in np.arange(0,total_number_of_planets,1):
							#### they'll be right next to each other!
							if fsk == 'ecs1':
								sampval = np.array(flat_samps['ecs'])[i][sample]						
							elif fsk == 'ecs2':
								sampval = np.array(flat_samps['ecs'])[i][sample]

							else:
								sampval = np.array(flat_samps[fsk])[i][sample] 
								#try:
								#	sampval = np.array(flat_samps[fsk])[i][sample]
								#except:
								#	try:
								#		sampval = np.array(flat_samps[fsk])[0][i][sample]		
								#	except:						
								#		sampval = np.array(flat_samps[fsk])[i][0][sample]

					#if (fsk != flat_samps_keys[-1]) and (i != total_number_of_planets - 1):
					if nfsk != len(flat_samps_keys) - 1:
						posterior_file.write(str(sampval)+',') #### there's still more to append

					elif nfsk == len(flat_samps_keys) - 1:
						posterior_file.write(str(sampval)+'\n') #### end of the line!
					else:
						print('SOMETHING WENT WRONG WRITING THE POSTERIOR FILE!')


				if continue_printing == 'y':
					print('sample ', sample)
					print('fsk = ', fsk)
					print('sampval = ', sampval)
					print(' ')


				#except:
				#	print('fsk = ', fsk)
				#	print(' ')
				#	traceback.print_exc()




				if continue_printing == 'y':
					continue_printing = input('Do you want to continue printing? y/n: ')



		posterior_file.close()
		print('Posteriors written to file.')



		if keep_showing == 'y':
			keep_showing = input('Do you want to keep showing plots? y/n: ')



	except:

		#error_log = open('/home/amteachey/Documents/Projects/Unistellar_transit_obs/kic_fitting_errors.txt', mode='a')
		error_log = open(moonpydir+'/fitting_errors.txt', mode='a')
		error_log.write('Exception raised for: '+str(self.target)+'\n')
		error_log.write('\n')
		error_log.close()

		traceback.print_exc()


		#continue_query = input('Exception was raised... continue? y/n: ')
		#if continue_query != 'y':
		#	raise Exception('you opted not to continue.')		

		#continue 
