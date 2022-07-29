from __future__ import division
import numpy as np
import astropy
from scipy.stats import norm,beta,truncnorm
import json
import matplotlib.pyplot as plt 
import traceback
import time
import os
import sys
from mp_batman import run_batman
import pandoramoon as pandora 
from pandoramoon.helpers import ld_convert, ld_invert 
from scipy.stats import norm 
try:
	import pyluna
except:
	print('Unable to load pyluna.')

try:
	from run_pandora import run_Pandora 
except:
	print('Unable to load run_pandora')

moonpydir = os.path.realpath(__file__)
moonpydir = moonpydir[:moonpydir.find('/mp_fit.py')]


### MULTINEST / ULTRANEST CUBE TRANSFORMATIONS
### note that in this context 'x' is the cube!
def transform_uniform(x,lower,upper):
	valrange = upper - lower 
	return lower + (valrange)*x

def transform_loguniform(x,lower,upper):
	l_lower=np.log10(lower)
	l_upper=np.log10(upper)
	valrange = l_upper - l_lower 
	#return np.exp(l_lower + x*(l_upper-l_lower))
	return 10**(l_lower + (valrange * x))


def transform_normal(x,mu,sigma):
	return norm.ppf(x,loc=mu,scale=sigma)

def transform_beta(x,lower,upper):
	return beta.ppf(x,lower,upper)

def transform_truncated_normal(x,mu,sigma,lower=0.,upper=1.):
	ar, br = (lower - mu) / sigma, (upper - mu) / sigma
	return truncnorm.ppf(x,ar,br,loc=mu,scale=sigma)

def transform_gauss(x, mu, sigma):
	#### redundant with transform_normal above
	gaussian = norm(loc=mu, scale=sigma)
	return gaussian.ppf(x)




### PYMULTINEST FUNCTIONS 
def pymn_prior(cube):
	#for pidx, parlabs, parprior, partuple in zip(np.arange(0,len(mn_prior_forms),1), mn_param_labels, mn_prior_forms, mn_limit_tuple):
	for pidx, parlabs, parprior, partuple in zip(np.arange(0,len(mn_variable_prior_forms),1), mn_variable_labels, mn_variable_prior_forms, mn_variable_limit_tuple):
		if parprior == 'uniform':
			cube[pidx] = transform_uniform(cube[pidx], partuple[0], partuple[1])
		elif parprior == 'loguniform':
			cube[pidx] = transform_loguniform(cube[pidx], partuple[0], partuple[1])
		elif parprior == 'normal':
			cube[pidx] = transform_normal(cube[pidx], partuple[0], partuple[1])
		elif parprior == 'beta':
			cube[pidx] = transform_beta(cube[pidx], partuple[0], partuple[1])
		elif parprior == 'truncnorm':
			cube[pidx] = transform_truncated_normal(cube[pidx], partuple[0], partuple[1], partuple[2], partuple[3])


def pymn_loglike_LUNA(cube):

	pymn_var_dict = {} #### dictionary of variables, will take the cube as the argument.
	pymn_fixed_dict = {} ### dictionary of fixed variables, will stay constant for each run.

	for pidx, parlab in enumerate(mn_variable_labels):
		pymn_var_dict[parlab] = cube[pidx]

	for pidx, parlab in enumerate(mn_fixed_labels):
		pymn_fixed_dict[parlab] = mn_param_dict[parlab][1] ### grabs the fixed value!

	### now you should be able to run_LUNA(param_dict)
	#LUNA_times, LUNA_fluxes = pyluna.run_LUNA(data_times, **pymn_param_dict, add_noise='n', show_plots='n')
	LUNA_times, LUNA_fluxes = pyluna.run_LUNA(data_times, **pymn_var_dict, **pymn_fixed_dict, model=mn_model, add_noise='n', show_plots='n')

	loglikelihood = np.nansum(-0.5 * ((LUNA_fluxes - data_fluxes) / data_errors)**2) ### SHOULD MAKE THIS BETTER, to super-penalize running out of bounds!
	return loglikelihood 


def pymn_loglike_batman(cube):

	pymn_var_dict = {} #### dictionary of variables, will take the cube as the argument.
	pymn_fixed_dict = {} ### dictionary of fixed variables, will stay constant for each run.

	for pidx, parlab in enumerate(mn_variable_labels):
		pymn_var_dict[parlab] = cube[pidx]

	for pidx, parlab in enumerate(mn_fixed_labels):
		pymn_fixed_dict[parlab] = mn_param_dict[parlab][1] ### grabs the fixed value!


	### now you should be able to run_LUNA(param_dict)
	batman_times, batman_fluxes = run_batman(data_times, **pymn_var_dict, **pymn_fixed_dict, model=mn_model, add_noise='n', show_plots='n')

	loglikelihood = np.nansum(-0.5 * ((batman_fluxes - data_fluxes) / data_errors)**2) ### SHOULD MAKE THIS BETTER, to super-penalize running out of bounds!
	return loglikelihood 



#### EMCEE FUNCTIONS
def emcee_lnprior(params):

	in_bounds = 'y'

	for param, parlab, parlim in zip(params, mc_variable_labels, mc_variable_limit_tuple):

		if (param < parlim[0]) or (param > parlim[1]):
			in_bounds = 'n'
			break
	if in_bounds == 'n':
		return -np.inf 
	else:
		return 0.0 


def emcee_lnlike_LUNA(params, data_times, data_fluxes, data_errors):
	emcee_param_dict = {} ### initialize it
	emcee_var_param_dict = {}
	emcee_fixed_param_dict = {}

	try:
		for pidx, parlab in enumerate(mc_variable_labels):
			emcee_var_param_dict[parlab] = params[pidx]

		for pidx, parlab in enumerate(mc_fixed_labels):
			emcee_fixed_param_dict[parlab] = mc_param_dict[parlab][1] ### fixed value!

		LUNA_times, LUNA_fluxes = pyluna.run_LUNA(data_times, **emcee_var_param_dict, **emcee_fixed_param_dict, model=mc_model, add_noise='n', show_plots='n')
	except:
		print('params = ', params)
		print(' ')
		print("emcee_var_param_dict: ", emcee_var_param_dict)
		print(' ')
		print('emcee_fixed_param_dict: ', emcee_fixed_param_dict)
		raise Exception('there was a critical failure of some kind.')


	loglikelihood = np.nansum(-0.5 * ((LUNA_fluxes - data_fluxes) / data_errors)**2) ### SHOULD MAKE THIS BETTER, to super-penalize running out of bounds!

	return loglikelihood 


def emcee_lnlike_batman(params, data_times, data_fluxes, data_errors):
	emcee_param_dict = {} ### initialize it
	emcee_var_param_dict = {}
	emcee_fixed_param_dict = {}

	for pidx, parlab in enumerate(mc_variable_labels):
		emcee_var_param_dict[parlab] = params[pidx]

	for pidx, parlab in enumerate(mc_fixed_labels):
		emcee_fixed_param_dict[parlab] = mc_param_dict[parlab][1] ### fixed value!

	batman_times, batman_fluxes = run_batman(data_times, **emcee_var_param_dict, **emcee_fixed_param_dict, model=mc_model, add_noise='n', show_plots='n')

	loglikelihood = np.nansum(-0.5 * ((batman_fluxes - data_fluxes) / data_errors)**2) ### SHOULD MAKE THIS BETTER, to super-penalize running out of bounds!

	return loglikelihood 


def emcee_lnprob_LUNA(params, data_times, data_fluxes, data_errors):
	lp = emcee_lnprior(params)
	if type(lp) != float:
		raise Exception('emcee_lnprior is not a float.')
	if not np.isfinite(lp):
		return -np.inf 
	return lp + emcee_lnlike_LUNA(params, data_times, data_fluxes, data_errors)


def emcee_lnprob_batman(params, data_times, data_fluxes, data_errors):
	lp = emcee_lnprior(params)
	if type(lp) != float:
		raise Exception('emcee_lnprior is not a float.')
	if not np.isfinite(lp):
		return -np.inf 
	return lp + emcee_lnlike_batman(params, data_times, data_fluxes, data_errors)







###### SAMPLERS BELOW !##########################


def mp_multinest(times, fluxes, errors, param_dict, nlive, targetID, model="M", nparams=14, modelcode='LUNA', show_plot='y'):
	import pymultinest
	### this function will start PyMultiNest!

	global mn_param_labels ### all inputs, fixed and variable
	global mn_fixed_labels ### all fixed inputs, not variable
	global mn_variable_labels ### al variable inputs, not fixed
	global mn_prior_forms
	global mn_variable_prior_forms
	global mn_fixed_prior_forms
	global mn_limit_tuple
	global mn_variable_limit_tuple
	global mn_fixed_limit_tuple
	global data_times
	global data_fluxes
	global data_errors
	global mn_param_dict
	global mn_model

	mn_param_labels = []
	mn_fixed_labels = []
	mn_variable_labels = []
	mn_prior_forms = []
	mn_fixed_prior_forms = []
	mn_variable_prior_forms = []
	mn_limit_tuple = []
	mn_fixed_limit_tuple = []
	mn_variable_limit_tuple = []
	data_times = []
	data_fluxes = []
	data_errors = []
	mn_param_dict = param_dict
	mn_model = model

	for parlab in param_dict.keys():
		parprior, partup = param_dict[parlab][0], param_dict[parlab][1]
		mn_param_labels.append(parlab)
		mn_prior_forms.append(parprior)
		mn_limit_tuple.append(partup)
		
		if parprior == 'fixed':
			mn_fixed_labels.append(parlab)
			mn_fixed_prior_forms.append(parprior)
			mn_fixed_limit_tuple.append(partup)
		else:
			mn_variable_labels.append(parlab)
			mn_variable_prior_forms.append(parprior)
			mn_variable_limit_tuple.append(partup)


	for t,f,e in zip(times, fluxes, errors):
		data_times.append(t)
		data_fluxes.append(f)
		data_errors.append(e)

	outputdir = moonpydir+'/MultiNest_fits'
	if os.path.exists(outputdir) == False:
		os.system('mkdir '+outputdir) ### create MultiNest_fits directory
	outputdir = outputdir+'/'+str(modelcode)
	if os.path.exists(outputdir) == False:
		os.system('mkdir '+outputdir) ### creates modelcode directory
	outputdir = outputdir+'/'+str(targetID)
	if os.path.exists(outputdir) == False:
		os.system('mkdir '+outputdir) ### creates target directory
	outputdir = outputdir+'/'+str(mn_model)
	if os.path.exists(outputdir) == False:
		os.system('mkdir '+outputdir) ### creates model directory
	outputdir = outputdir+'/chains'
	if os.path.exists(outputdir) == False:
		os.system('mkdir '+outputdir) ### creates chains directory


	if modelcode == 'LUNA':
		pymultinest.run(LogLikelihood=pymn_loglike_LUNA, Prior=pymn_prior, n_dims=nparams, n_live_points=nlive, outputfiles_basename=outputdir+'/'+str(targetID), resume=True, verbose=True)
	elif modelcode == "batman":
		pymultinest.run(LogLikelihood=pymn_loglike_batman, Prior=pymn_prior, n_dims=nparams, n_live_points=nlive, outputfiles_basename=outputdir+'/'+str(targetID), resume=True, verbose=True)
	
	json.dump(param_labels, open(outputdir+'/'+str(targetID)+"_params.json", 'w')) ### save parameter names






def mp_ultranest(times, fluxes, errors, param_dict, nlive, targetID, model="M", modelcode='Pandora', resume='resume', show_plot='y'):	
	import ultranest 
	import ultranest.stepsampler 
	from ultranest import ReactiveNestedSampler 
	from ultranest.plot import cornerplot 


	### this function will start UltraNest 


	#### DEFINE THESE GLOBALS, BECAUSE THEY WILL BE NEEDED -- I FORGOT WHY, BUT TRUST ME.
	global un_param_labels ### all inputs, fixed and variable
	global un_fixed_labels ### all fixed inputs, not variable
	global un_variable_labels ### al variable inputs, not fixed
	global un_prior_forms
	global un_variable_prior_forms
	global un_fixed_prior_forms
	global un_limit_tuple
	global un_variable_limit_tuple
	global un_fixed_limit_tuple
	global data_times
	global data_fluxes
	global data_errors
	global un_param_dict
	global un_model
	global nepochs 

	### labels
	un_param_labels = []
	un_fixed_labels = []
	un_variable_labels = []

	#### prior forms
	un_prior_forms = []
	un_fixed_prior_forms = []
	un_variable_prior_forms = []
	wrapped_params = []

	#### limit tuples
	un_limit_tuple = []
	un_fixed_limit_tuple = []
	un_variable_limit_tuple = []

	data_times = []
	data_fluxes = []
	data_errors = []

	un_param_dict = param_dict
	un_model = model

	#for parlab, parprior, partup in zip(param_labels, param_prior_forms, param_limit_tuple):
	for parlab in param_dict.keys():
		parprior, partup = param_dict[parlab][0], param_dict[parlab][1]
		un_param_labels.append(parlab)
		un_prior_forms.append(parprior)
		un_limit_tuple.append(partup)
		
		if parprior == 'fixed':
			un_fixed_labels.append(parlab)
			un_fixed_prior_forms.append(parprior)
			un_fixed_limit_tuple.append(partup)

		else:
			un_variable_labels.append(parlab)
			un_variable_prior_forms.append(parprior)
			un_variable_limit_tuple.append(partup)
		
		#if parlab == 'tau':
		#	wrapped_params.append(True)
		#else:
		#	wrapped_params.append(False)



	for t,f,e in zip(times, fluxes, errors):
		data_times.append(t)
		data_fluxes.append(f)
		data_errors.append(e)
	data_times, data_fluxes, data_errors = np.array(data_times), np.array(data_fluxes), np.array(data_errors)


	print('param_dict.keys() = ', param_dict.keys())

	#### reduce the times to just the ones in the vicinty of the transits.
	all_taus = []
	if param_dict['t0_bary'][0] == 'fixed':
		first_tau = float(param_dict['t0_bary'][1])
	else:
		first_tau = float(param_dict['t0_bary'][1][0]) ### list includes the prior type, the expected answer and the sigma in all cases
	per_bary = float(param_dict['per_bary'][1][0])


	while np.nanmin(data_times) + per_bary < first_tau: ### implies first_tau is more than one period from the start of the baseline
		first_tau = first_tau - per_bary 
	all_taus.append(first_tau)

	while all_taus[-1] < np.nanmax(data_times):
		all_taus.append(all_taus[-1] + per_bary)

	nepochs = 0
	epoch_duration = 5 #### days

	for ntau, tau in enumerate(all_taus):
		#### grab the indices within the desired distance.
		tau_idxs = np.where((data_times > tau - 0.5*epoch_duration) & (data_times <= tau + 0.5*epoch_duration))[0]
		tau_times = data_times[tau_idxs]

		if len(tau_times) > 0:
			#### this epoch is in the data
			nepochs += 1

		if ntau == 0:
			model_idxs = tau_idxs 
			model_times = tau_times
		else:
			model_idxs = np.concatenate((model_idxs, tau_idxs))
			model_times = np.concatenate((model_times, tau_times))	



	#### THESE ARE THE DATA YOU ARE GOING TO MODEL!!!
	data_times, data_fluxes, data_errors = np.array(data_times)[model_idxs], np.array(data_fluxes)[model_idxs], np.array(data_errors)[model_idxs]

	plt.scatter(data_times, data_fluxes, s=20, facecolor='DodgerBlue', edgecolor='k')
	plt.title('Data to be fit with '+str(modelcode)+' and UltraNest')
	plt.show()

	print('NUMBER OF DATA POINTS TO FIT: ', len(data_times))

	### MAKE THE OUTPUT DIRECTORY IF IT DOESN'T ALREADY EXIST
	outputdir = moonpydir+'/UltraNest_fits'
	if os.path.exists(outputdir) == False:
		os.system('mkdir '+outputdir) ### create UltraNest_fits directory

	outputdir = outputdir+'/'+str(modelcode)
	if os.path.exists(outputdir) == False:
		os.system('mkdir '+outputdir) ### creates modelcode directory

	outputdir = outputdir+'/'+str(targetID)
	if os.path.exists(outputdir) == False:
		os.system('mkdir '+outputdir) ### creates target directory

	outputdir = outputdir+'/'+str(model)
	print('model: ', model)
	if os.path.exists(outputdir) == False:
		os.system('mkdir '+outputdir) ### creates model directory

	print('outputdir: ', outputdir)


	if modelcode == 'Pandora':

		#### DEFINE HERE FOR THE SAKE OF AVOIDING DICTIONARY CALLS -- ALTHOUGH THIS DOESN'T SEEM TO BE SAVING MUCH TIME.
		prior_mins, prior_maxes, prior_mus, prior_sigmas, prior_types = [], [], [], [], []

		#prior_ranges, prior_mins = [], []
		for pidx, parlabs, parprior, partuple in zip(np.arange(0,len(un_variable_prior_forms),1), un_variable_labels, un_variable_prior_forms, un_variable_limit_tuple):
			if (parprior == 'uniform') or (parprior == 'loguniform'):
				prior_mins.append(partuple[0])
				prior_maxes.append(partuple[1])
				prior_mus.append(None)
				prior_sigmas.append(None)
				prior_types.append(parprior)

			elif (parprior == 'normal') or (parprior == 'gaussian'):
				prior_mins.append(None)
				prior_maxes.append(None)
				prior_mus.append(partuple[0])
				prior_sigmas.append(partuple[1])
				prior_types.append(parprior)

			elif parprior == 'truncnormal':
				prior_mins.append(partuple[2])
				prior_maxes.append(partuple[3])
				prior_mus.append(partuple[0])
				prior_sigmas.append(partuple[1])
				prior_types.append(parprior)

		prior_mins, prior_maxes, prior_mus, prior_sigmas, prior_types = np.array(prior_mins), np.array(prior_maxes), np.array(prior_mus), np.array(prior_sigmas), np.array(prior_types)


		def ultn_transform(cube):
			p = cube.copy()

			for pidx in np.arange(0,len(prior_types),1):
				if prior_types[pidx] == 'uniform':
					p[pidx] = transform_uniform(x=cube[pidx], lower=prior_mins[pidx], upper=prior_maxes[pidx])
				elif prior_types[pidx] == 'loguniform':
					p[pidx] = transform_loguniform(x=cube[pidx], lower=prior_mins[pidx], upper=prior_maxes[pidx])
				elif (prior_types[pidx] == 'normal') or (prior_types[pidx] == 'gaussian'):
					p[pidx] = transform_normal(x=cube[pidx], mu=prior_mus[pidx], sigma=prior_sigmas[pidx]) 
				elif (prior_types[pidx] == 'truncnorm'):
					p[pidx] = transform_truncated_normal(x=cube[pidx], mu=prior_mus[pidx], sigma=prior_sigmas[pidx], lower=prior_mins[pidx], upper=prior_maxes[pidx])

			return p

		def ultn_loglike_Pandora(cube):

			#### THIS FUNCTION GETS CALLED REPEATEDLY! SO
			##### YOU WALK AROUND IN q1,q2 but model in u1,u2 at every step.

			ultn_var_dict = {} #### dictionary of variables, will take the cube as the argument.
			ultn_fixed_dict = {} ### dictionary of fixed variables, will stay constant for each run.

			for pidx, parlab in enumerate(un_variable_labels):
				ultn_var_dict[parlab] = cube[pidx]

			for pidx, parlab in enumerate(un_fixed_labels):
				ultn_fixed_dict[parlab] = un_param_dict[parlab][1] ### grabs the fixed value!

			q1, q2 = ultn_var_dict['q1'], ultn_var_dict['q2']
			u1, u2 = ld_convert(q1, q2)
			
			#### add these to the dictionary
			ultn_var_dict['u1'] = u1
			ultn_var_dict['u2'] = u2 
			#### remove these from the dictionary 
			ultn_var_dict.pop('q1', None)
			ultn_var_dict.pop('q2', None) 
			

			flux_trial_planet, flux_trial_moon, flux_trial_total, xp_trial, yp_trial, xm_trial, ym_trial = pandora.pandora(**ultn_var_dict, **ultn_fixed_dict, 
		   		epoch_distance = ultn_var_dict['per_bary'],
		   		supersampling_factor = 1, 
		   		occult_small_threshold = 0.01, 
		   		hill_sphere_threshold=1.0,
		   		numerical_grid=25,
		   		time=data_times, 
		   		)

			loglike = -0.5 * np.nansum(((flux_trial_total - data_fluxes) / data_errors)**2)

			return loglike


		#### PANDORA IS CALLED by calling ultn_loglike_Pandora below
		try:
			#sampler = ReactiveNestedSampler(un_variable_labels, ultn_loglike_Pandora, transform=ultn_transform, log_dir=outputdir, resume='resume', vectorized=True)	
			sampler = ReactiveNestedSampler(un_variable_labels, ultn_loglike_Pandora, transform=ultn_transform, log_dir=outputdir, resume='resume')										
			sampler.stepsampler = ultranest.stepsampler.RegionSliceSampler(nsteps=4000, adaptive_nsteps='move-distance')
			result_planet_moon = sampler.run(min_num_live_points=nlive)


		except AssertionError:
			#try:
			os.system('rm -rf '+outputdir)
			os.system('mkdir '+outputdir) ### creates model directory
			print('COULD NOT RUN THE SAMPLER WITH RESUME=RESUME. TRYING RESUME=OVERWRITE.')
			sampler = ReactiveNestedSampler(un_variable_labels, ultn_loglike_Pandora, transform=ultn_transform, log_dir=outputdir, resume='overwrite')	
			sampler.stepsampler = ultranest.stepsampler.RegionSliceSampler(nsteps=4000, adaptive_nsteps='move-distance')		
			result_planet_moon = sampler.run(min_num_live_points=nlive)
		
		try:
			sampler.print_results()
		except:
			traceback.print_exc()
			print(' ')
			print("WARNING: could not call sampler.print_results().")
			time.sleep(3)
		

		cornerplot(result_planet_moon)
		plt.savefig(outputdir+'/model'+str(model)+'_cornerplot.png', dpi=300)
		plt.show()

		try:
			sampler.plot()
		except:
			traceback.print_exc()
			print(' ')
			print('WARNING: could not call sampler.plot().')
			time.sleep(3)

		pew = np.genfromtxt(outputdir+'/chains/equal_weighted_post.txt')
		pewtxt = open(outputdir+'/chains/equal_weighted_post.txt', mode='r')

		firstline = pewtxt.readline()
		cols = firstline.split()
		pewtxt.close()
		pewdict = {}
		for ncol,col in enumerate(cols):
			pewdict[col] = pew.T[ncol][1:]
			n,bins,edges = plt.hist(pewdict[col], bins=20, facecolor='DodgerBlue', edgecolor='k')
			plt.xlabel(col)
			plt.show()

		sampler.plot_trace()

	json.dump(un_param_labels, open(outputdir+'/'+str(targetID)+"_model"+str(model)+"_params.json", 'w')) ### save parameter names








def mp_emcee(times, fluxes, errors, param_dict, nwalkers, nsteps, targetID, nparams=14, resume=True, modelcode="LUNA", model='M', storechain=False, burnin_pct=0.1, show_plot='y'):
	import emcee
	import corner

	global mc_param_labels ### all inputs
	global mc_fixed_labels ### all fixed inputs, not variable
	global mc_variable_labels ### al variable inputs, not fixex
	global mc_prior_forms
	global mc_variable_prior_forms
	global mc_fixed_prior_forms
	global mc_limit_tuple
	global mc_variable_limit_tuple
	global mc_fixed_limit_tuple
	global data_times
	global data_fluxes
	global data_errors
	global mc_param_dict
	global mc_model

	mc_param_labels = []
	mc_fixed_labels = []
	mc_variable_labels = []
	mc_prior_forms = []
	mc_fixed_prior_forms = []
	mc_variable_prior_forms = []
	mc_limit_tuple = []
	mc_fixed_limit_tuple = []
	mc_variable_limit_tuple = []
	data_times = []
	data_fluxes = []
	data_errors = []
	mc_param_dict = param_dict
	mc_model = model

	#for parlab, parprior, partup in zip(param_labels, param_prior_forms, param_limit_tuple):
	for parlab in param_dict.keys():
		parprior, partup = param_dict[parlab][0], param_dict[parlab][1]
		mc_param_labels.append(parlab)
		mc_prior_forms.append(parprior)
		mc_limit_tuple.append(partup)
		
		if parprior == 'fixed':
			mc_fixed_labels.append(parlab)
			mc_fixed_prior_forms.append(parprior)
			mc_fixed_limit_tuple.append(partup)
		else:
			mc_variable_labels.append(parlab)
			mc_variable_prior_forms.append(parprior)
			mc_variable_limit_tuple.append(partup)

	for t,f,e in zip(times, fluxes, errors):
		data_times.append(t)
		data_fluxes.append(f)
		data_errors.append(e)


	outputdir = 'emcee_fits'
	if os.path.exists(outputdir) == False:
		os.system('mkdir '+outputdir)
	outputdir = outputdir+'/'+str(modelcode)
	if os.path.exists(outputdir) == False:
		os.system('mkdir '+outputdir)
	outputdir = outputdir+'/'+str(targetID)
	if os.path.exists(outputdir) == False:
		os.system('mkdir '+outputdir)
	outputdir = outputdir+'/'+str(mc_model)
	if os.path.exists(outputdir) == False:
		os.system('mkdir '+outputdir)
	outputdir = outputdir+'/chains'
	if os.path.exists(outputdir) == False:
		os.system('mkdir '+outputdir)


	chainfile = outputdir+'/'+str(targetID)+'_mcmc_chain.txt'
	lnpostfile = outputdir+'/'+str(targetID)+"_mcmc_lnpost.txt"


	### if chainfile doesn't exist, resume == False!
	if os.path.exists(chainfile) == False:
		resume = False 

	### INITIALIZE EMCEE PARAMETERS
	ndim = len(mc_variable_labels)
	assert ndim == nparams
	#pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
	### pos is a list of arrays!
	if resume == False:
		chainf = open(chainfile, mode='w')
		chainf.close()

		lnpostf = open(lnpostfile, mode='w')
		lnpostf.close()

		mcmc_iter = nsteps 

		### initialize the positions!
		pos = []
		for walker in np.arange(0,nwalkers,1):
			walker_pos = []
			for partup in mc_variable_limit_tuple:
				if np.abs(partup[1] - partup[0]) > 1e3: ### implies huge jumps:
					parspot = np.random.choice(np.logspace(start=np.log10(partup[0]), stop=np.log10(partup[1]), num=1e4))
				else:
					parspot = np.random.choice(np.linspace(partup[0], partup[1], 1e4))

				walker_pos.append(parspot)
			walker_pos = np.array(walker_pos)
			pos.append(walker_pos)

		#print('pos = ', pos)


	elif resume == True:
		### load the positions!
		paused = np.genfromtxt(chainfile)
		pos = []
		### the last N=nwalker lines should be the positions of the walkers!
		for walkerline in paused[-nwalkers:]:
			pos.append(np.array(walkerline[1:]))
		#pos = paused[-nwalkers, -ndim:]

		#print('pos = ', pos)
		### MAY 29 -- the problem with resume is that pos does not have the sahem shape as pos above!
		### when resume == False, pos is a list of arrays... presumably 100 arrays (= nwalkers), each with len=ndim.
		### but in the above sample, pos is only a single array with len=ndim.
		### needs to be identical to the above (which works!) in order to work correctly.

		mcmc_iter = int((nsteps - (np.shape(paused)[0]/nwalkers)))


	#print('pos = ', pos)
	#print(' ')
	#continue_query = input('Continue? y/n: ')
	#if continue_query != 'y':
	#	raise Exception('you opted not to continue.')

	if modelcode == 'LUNA':
		print('running the sampler with LUNA.')
		print('nwalkers = ', nwalkers)
		print("ndim = ", ndim)
		sampler = emcee.EnsembleSampler(nwalkers, ndim, emcee_lnprob_LUNA, args=(data_times, data_fluxes, data_errors))

	elif modelcode == "batman":
		print('running the sampler with batman.')
		sampler = emcee.EnsembleSampler(nwalkers, ndim, emcee_lnprob_batman, args=(data_times, data_fluxes, data_errors))

	### run the sampler 
	#sampler.run_mcmc(pos, 10000, lnprob0=np.linspace(0,0,nwalkers))

	#sampler.run_mcmc(pos, 10000)
	for residx, result in enumerate(sampler.sample(pos, iterations=mcmc_iter, storechain=storechain)):

		### WRITE OUT CHAIN TO A FILE.
		position, lnprob, rstate = result[0], result[1], result[2]

		### write the walker position ot the chainfile
		chainf = open(chainfile, mode='a')
		for k in range(position.shape[0]):
			chainf.write("{0:4d} {1:s}\n".format(k, " ".join(str(q) for q in position[k])))
		chainf.close()

		lnpostf = open(lnpostfile, mode='a')
		for i in range(lnprob.shape[0]):
			lnpostf.write("{0:4d} {1:s}\n".format(i, str(lnprob[i])))
		lnpostf.close()

		### progress bar:
		print('step '+str(residx)+' of '+str(mcmc_iter))
		print(' ')
		#width=30
		#n = int((width+1) * float(residx) / mcmc_iter)
		#sys.stdout.write("\r[{0}{1}]".format('#' * n, ' ' * (width - n)))





	#for i, result in enumerate(sampler.sample(p0, iterations=10000)):
	#	if (i+1) % 100 == 0:
	#	   print("{0:5.1%}".format(float(i) / nsteps))
	if storechain == True:
		### note that 'samples' should have shape = ((nsteps - n_burnin) * nwalkers, ndim)
		### example: if you have nwalkers=30, n_burning=10, nsteps=200, and ndim=11,
		### your final shape should be (5700, 11), or (200 - 10) * nwalkers, ndim.
		print('sampler.chain.shape = ', sampler.chain.shape)
		#samples = sampler.chain[:,10:,:].reshape((-1, ndim))
		samples = sampler.chain[:,burnin_pct*mcmc_iter,:].reshape((-1,ndim))
		print('sampler.chain.reshaped.shape = ', sampler.chain[:,10:,:].reshape((-1,ndim)).shape)
	else:
		### load the chainsfile!
		samples = np.genfromtxt(chainfile)[(nsteps - (int(burnin_pct*nsteps))):,1:]

	if show_plot == 'y':
		fig = corner.corner(samples, labels=mc_variable_labels)
		plt.savefig(outputdir+'/'+str(targetID)+"_corner.png")
		plt.close()


