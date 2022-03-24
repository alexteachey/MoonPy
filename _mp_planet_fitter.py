from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pandas
import os
import traceback
from astropy.io import fits
import exoplanet as xo
import astropy.units as u 
import pymc3 as pm
import pymc3_ext as pmx
import arviz as az
try:
	import corner
except:
	print('could not import corner.')
from celerite2.theano import terms, GaussianProcess
import aesara_theano_fallback.tensor as tt
import platform 
from matplotlib import rcParams
from moonpy import * 

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


def run_planet_fit(self, period=None, tau0=None, tdur_hours=None, smass=None, show_plots=True, use_mp_detrend=False):

	if show_plots == True:
		keep_showing = 'y'
	else:
		keep_showing = 'n'

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
	except:
		manual_smass_entry = input('Enter stellar mass, in solar units: ')
		self.smass = float(manual_smass_entry)
		






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
		near_transit_idxs = []
		for ncctime, cctime in enumerate(cctimes):
			if np.any(np.abs(cctime - taus) < 10):
				near_transit_idxs.append(ncctime)
		near_transit_idxs = np.array(near_transit_idxs)

		print('len(cctimes) = ', len(cctimes))
		print('len(near_transit_idxs) = ', len(near_transit_idxs))
		continue_query = input('Do you wish to continue? y/n: ')
		if continue_query != 'y':
			raise Exception('you opted not to continue.')

		cctimes, ccsap, ccsap_err, cc_fluxes_detrend, cc_errors_detrend, ccflags = cctimes[near_transit_idxs], ccsap[near_transit_idxs], ccsap_err[near_transit_idxs], cc_fluxes_detrend[near_transit_idxs], cc_errors_detrend[near_transit_idxs], ccflags[near_transit_idxs] 
		### update this
		good_flags = np.where(ccflags == 0)[0]


		#### test plot 
		#plt.scatter(cctimes, cc_fluxes_detrend, marker='o', facecolor='LightCoral', edgecolor='k', s=20)
		#plt.show()



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

			ax[0].scatter(cctimes[good_flags], cc_fluxes_detrend[good_flags], facecolor='LightCoral', edgecolor='k', s=20, alpha=0.7)

			for btt in transit_times:
				#### mark with a vertical line
				ax[0].plot(np.linspace(btt,btt,100), np.linspace(0.95*np.nanmin(cc_fluxes_detrend), 1.05*np.nanmax(cc_fluxes_detrend), 100), color='red', alpha=0.7, linewidth=1)
				#plt.plot(np.linspace(btt-0.5*duration,btt-0.5*duration,100), np.linspace(0.95*np.nanmin(cc_fluxes_detrend), 1.05*np.nanmax(cc_fluxes_detrend), 100), color='red', alpha=0.8, linestyle='--')
				#plt.plot(np.linspace(btt+0.5*duration,btt+0.5*duration,100), np.linspace(0.95*np.nanmin(cc_fluxes_detrend), 1.05*np.nanmax(cc_fluxes_detrend), 100), color='red', alpha=0.8, linestyle='--')			

			ax[1].scatter(fold_times[good_flags], ccksp[good_flags], facecolor='LightCoral', edgecolor='k', s=20, alpha=0.7)

			ax[0].set_title(kic)
			#ax[1].set_xlabel('BTJD')
			ax[0].set_ylabel('FLUX')
			ax[1].set_ylabel('FLUX')
			plt.show()




		#### NOW THE TRANSIT FITTING!
		#### The transit model in PyMC3

		#### MY TRY
		periods = self.period #### my guess
		t0s = first_transit
		model_times = cctimes[good_flags]
		#yerr = np.nanmedian(np.concatenate(kspsap_errors)[good_flags])
		#yerr = np.nanmedian(np.concatenate(errors_detrend)[good_flags])

		#if use_mp_detrend == True:
		yvals = cc_fluxes_detrend[good_flags]
		yerr = cc_errors_detrend[good_flags]
		
		#elif use_mp_detrend == False:
		#	#### normalize them 
		#	yvals = ccsap[good_flags] / np.nanmedian(ccsap[good_flags])
		#	yerr = ccsap_err[good_flags] / ccsap[good_flags]


		"""
		YOU SHOULD REALLY FOLLOW THIS PAGE -- TO INCORPORATE THE GPs!!!
		https://gallery.exoplanet.codes/tutorials/tess/#the-transit-model-in-pymc3

		"""


		with pm.Model() as model: #### this is a PYMC3 MODEL! 


			phase_lc = np.linspace(-0.3, 0.3, 100) #### what is this?!

			#### let y be the real data
			#yvals = np.concatenate(kspsap_fluxes)[good_flags]
			#yvals = cc_fluxes_detrend[good_flags] 

			#### COLLECTING ALL THE PRIORS HERE.
			#mean = pm.Normal("mean", mu=1.0, sd=np.nanstd(np.concatenate(fluxes)[good_flags]))		
			mean = pm.Normal("mean", mu=1.0, sd=np.nanmedian(cc_errors_detrend))
			#ldcs = xo.distributions.QuadLimbDark("ldcs", testval=np.array([0.3, 0.2])) #### STICK WITH THIS I GUESS		
			ldcs = xo.distributions.QuadLimbDark("ldcs", testval=(0.3, 0.2)) #### making it a tuple -- I guess this is what they want?
			star = xo.LimbDarkLightCurve(ldcs[0], ldcs[1])

			BoundedNormal = pm.Bound(pm.Normal, lower=0.0)
			#m_star = BoundedNormal("m_star", mu=QLP_MASS, sd=0.05*QLP_MASS) ### now uncertainties available, just use 5%
			try:
				m_star = BoundedNormal("m_star", mu=self.smass, sd=np.nanmean(np.abs(self.smass_err)))
			except:
				m_star = BoundedNormal("m_star", mu=self.smass, sd=0.05*self.smass)
			#r_star = BoundedNormal("r_star", mu=QLP_RADIUS, sd=0.05*QLP_RADIUS) ### no uncertainties available, just using 5%
			r_star = BoundedNormal("r_star", mu=self.rstar_rsol, sd=0.05*self.rstar_rsol)


			t0 = pm.Normal("t0", mu=t0s, sd=1.0, shape=1)
			log10Period = pm.Normal("log10Period", mu=np.log10(periods), sd=0.1, shape=1)
			period = pm.Deterministic("period", 10**log10Period)	

			impact = pm.Uniform("impact", lower=0, upper=1) 
			log_depth = pm.Normal("log_depth", mu=np.log(np.nanmin(yvals)), sigma=2.0) ### why sigma=2?

			#rprstar = pm.Uniform("rprstar", lower=0, upper=1, shape=1, testval=np.array([0.05])) ### MY VERSION
			try:
				rprstar = pm.Uniform('rprstar', lower=0, upper=1, shape=1, testval=self.rprstar)
			except:
				#rprstar = pm.Uniform("rprstar", lower=0, upper=1, shape=1, testval=np.array([0.05])) ### MY VERSION	
				rprstar = pm.Uniform('rprstar', lower=0, upper=1, shape=1, testval=0.05)

			rplan = pm.Deterministic("rplan", rprstar * r_star) 
			#ecs = pmx.UnitDisk('ecs', testval=np.array([0.01, 0.0])) #### WHAT IS THIS?? 
			ecs = pmx.UnitDisk('ecs', testval=(0.01, 0.0)) #### making it a tuple... I think this is what they want.
			#ecc = xo.distributions.eccentricity.kipping13("ecc", lower=None, upper=None) ### NEW	
			ecc = pm.Deterministic("ecc", tt.sum(ecs ** 2))					
			#omega = pm.Uniform('omega', lower=0, upper=2*np.pi, shape=1, testval=np.array([np.pi])) ### NEW
			omega = pm.Deterministic("omega", tt.arctan2(ecs[1], ecs[0])) #### what is this??!
			xo.eccentricity.kipping13("ecc_prior", fixed=True, observed=ecc)
			
			#### GET RID OF THIS IF YOU'RE GOING TO USE THE QLP MASS AND RADIUS VALUES!
			#rhostar = pm.HalfNormal('rhostar', sigma=5) ### g/cm^3 ### based on the distribution on NASA Exoplanet Archive (eyeballed)

			#Transit jitter & GP parameters
			log_sigma_lc = pm.Normal("log_sigma_lc", mu=np.log(np.nanstd(yvals)), sd=10) #### why sd=10?
			log_rho_gp = pm.Normal("log_rho_gp", mu=0, sd=10) #### what is this mean and sd?!?!
			log_sigma_gp = pm.Normal("log_sigma_gp", mu=np.log(np.nanstd(yvals)), sd=10)  #### 


			# Set up a Keplerian orbit for the planets -- models the system
			orbit = xo.orbits.KeplerianOrbit(period=period, t0=t0, b=impact, ecc=ecc, omega=omega, m_star=m_star, r_star=r_star)

			#### creates an observation from the model
			light_curves = star.get_light_curve(orbit=orbit, r=rplan, t=model_times)	
			light_curve = tt.sum(light_curves, axis=-1) + mean 

			residuals = yvals - light_curve 


			### GP MODEL FOR THE LIGHT CURVE
			kernel = terms.SHOTerm(sigma=tt.exp(log_sigma_gp), rho=tt.exp(log_rho_gp), Q=1 / np.sqrt(2),)
			gp = GaussianProcess(kernel, t=model_times, yerr=tt.exp(log_sigma_lc))
			gp.marginal("gp", observed=residuals)
			



			# Here we track the value of the model light curve for plotting
			# purposes
			pm.Deterministic("light_curves", light_curves)
			pm.Deterministic("lc_pred", 1e6 * star.get_light_curve(orbit=orbit, r=rplan, t=t0 + phase_lc)[...,0]) #### WHAT IS ALL THIS?!



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
			#map_soln = pmx.optimize(start=map_soln, vars=(log10Period, t0))
			map_soln = pmx.optimize(start=map_soln, vars=[ldcs])
			map_soln = pmx.optimize(start=map_soln, vars=[log_depth]) ### huh? why again?
			map_soln = pmx.optimize(start=map_soln, vars=[impact]) #### why again?
			map_soln = pmx.optimize(start=map_soln, vars=[ecs])
			map_soln = pmx.optimize(start=map_soln, vars=[mean])
			map_soln = pmx.optimize(start=map_soln, vars=[log_sigma_lc, log_sigma_gp, log_rho_gp]) #huh?
			map_soln = pmx.optimize(start=map_soln) 
			

			#### RE-DOING WITH TUPLES
			"""
			map_soln = pmx.optimize(start=model.test_point, vars=np.array([(log_sigma_lc, log_sigma_gp, log_rho_gp)]))
			map_soln = pmx.optimize(start=map_soln, vars=np.array([(log_depth)]))
			map_soln = pmx.optimize(start=map_soln, vars=np.array([(impact)]))
			map_soln = pmx.optimize(start=map_soln, vars=np.array([(log10Period, t0)]))
			map_soln = pmx.optimize(start=map_soln, vars=np.array([(log10Period, t0)]))
			map_soln = pmx.optimize(start=map_soln, vars=np.array([(ldcs)]))
			map_soln = pmx.optimize(start=map_soln, vars=np.array([(log_depth)])) ### huh? why again?
			map_soln = pmx.optimize(start=map_soln, vars=np.array([(impact)])) #### why again?
			map_soln = pmx.optimize(start=map_soln, vars=np.array([(ecs)]))
			map_soln = pmx.optimize(start=map_soln, vars=np.array([(mean)]))
			map_soln = pmx.optimize(start=map_soln, vars=np.array([(log_sigma_lc, log_sigma_gp, log_rho_gp)])) #huh?
			map_soln = pmx.optimize(start=map_soln) 
			"""


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
		ax.scatter(model_times, (yvals - gp_mod)*1e6, facecolor='LightCoral', edgecolor='k', marker='o', s=20, alpha=0.5, label="de-trended data", zorder=0)
		mod = extras["light_curves"][:,0]
		scatter_ppm = np.nanstd((yvals - mod)*1e6)

		ax.plot(model_times, mod*1e6, label='transit model', zorder=1)
		ax.set_ylim(np.nanmin((yvals-gp_mod)*1e6) - scatter_ppm, 2*scatter_ppm)
		ax.legend(fontsize=10, loc='upper right')
		ax.set_ylabel("de-trended flux [ppm]")


		ax = axes[2]
		mod = gp_mod + np.sum(extras["light_curves"], axis=-1)
		ax.scatter(model_times, (yvals - mod)*1e6, facecolor='LightCoral', edgecolor='k', marker='o', s=20, alpha=0.5, zorder=0)
		ax.axhline(0, color="#aaaaaa", lw=1)
		ax.set_ylabel("residuals [ppm]")
		ax.set_xlim(model_times.min(), model_times.max())
		ax.set_ylim(-2*scatter_ppm, 2*scatter_ppm)
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

		# Compute the GP prediction
		gp_mod = extras["gp_pred"] + map_soln["mean"]  # np.median(
		#     flat_samps["gp_pred"].values + flat_samps["mean"].values[None, :], axis=-1
		# )

		# Get the posterior median orbital parameters
		post_period = np.median(flat_samps["period"])
		post_t0 = np.median(flat_samps["t0"])

		# Plot the folded data
		x_fold = (model_times - post_t0 + 0.5 * post_period) % post_period - 0.5 * post_period
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
		plt.xlim(-0.5 * post_period, 0.5 * post_period)
		plt.xlabel("time since transit [days]")
		plt.ylabel("de-trended flux [ppm]")
		_ = plt.xlim(-0.15, 0.15)
		_ = plt.ylim(np.nanmin(yvals_ppm)-scatter_ppm, scatter_ppm)

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




		flat_samps_keys = ['mean','ldc1','ldc2','m_star','r_star','t0','log10Period','period','impact','log_depth','rprstar','rplan','ecs1','ecs2','ecc','omega','log_sigma_lc','log_rho_gp','log_sigma_gp']
		#posterior_file = open(kicpath+'/'+kic+'_posteriors.csv', mode='w')
		posterior_file = open(savepath+'/'+self.target+'_posteriors.csv', mode='w')
		
		### write the header
		for nfsk,fsk in enumerate(flat_samps_keys):
			if fsk != flat_samps_keys[-1]:
				posterior_file.write(fsk+',')
			elif fsk == flat_samps_keys[-1]:
				posterior_file.write(fsk+'\n')

		extra_idx_keys = ['t0', 'log10Period', 'period', 'rprstar', 'rplan']

		for sample in np.arange(0,nsamples,1):
			for fsk in flat_samps_keys:
				if fsk == 'ldc1':
					sampval = np.array(flat_samps['ldcs'])[0][sample]
				elif fsk == 'ldc2':
					sampval = np.array(flat_samps['ldcs'])[1][sample]	
				elif fsk == 'ecs1':
					sampval = np.array(flat_samps['ecs'])[0][sample]						
				elif fsk == 'ecs2':
					sampval = np.array(flat_samps['ecs'])[1][sample]	

				elif fsk in extra_idx_keys:
					sampval = np.array(flat_samps[fsk])[0][sample] #### have to do this for some stupid reason.

				else:
					sampval = np.array(flat_samps[fsk])[sample]
				
				if fsk != flat_samps_keys[-1]:
					posterior_file.write(str(sampval)+',')
				elif fsk == flat_samps_keys[-1]:
					posterior_file.write(str(sampval)+'\n')

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