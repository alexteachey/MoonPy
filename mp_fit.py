from __future__ import division
import numpy as np
import astropy
from scipy.stats import norm,beta,truncnorm
import pyluna
import pymultinest
import json
import matplotlib.pyplot as plt 
import os

### MULTINEST CUBE TRANSFORMATIONS
### note that in this context 'x' is the cube!
def transform_uniform(x,a,b):
    return a + (b-a)*x

def transform_loguniform(x,a,b):
    la=np.log(a)
    lb=np.log(b)
    return np.exp(la + x*(lb-la))

def transform_normal(x,mu,sigma):
    return norm.ppf(x,loc=mu,scale=sigma)

def transform_beta(x,a,b):
    return beta.ppf(x,a,b)

def transform_truncated_normal(x,mu,sigma,a=0.,b=1.):
    ar, br = (a - mu) / sigma, (b - mu) / sigma
    return truncnorm.ppf(x,ar,br,loc=mu,scale=sigma)


def pymn_prior(cube, ndim, nparams):
	### hopefully you can inherit these within the mp_multinest function!
	#for pidx, parprior, partuple in zip(np.arange(0,len(param_prior_forms),1), param_prior_forms, param_limit_tuple):
	for pidx, parprior, partuple in zip(np.arange(0,len(mn_prior_forms),1), mn_prior_forms, mn_limit_tuple):
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

def pymn_loglike_LUNA(cube, ndim, nparams):
	### have to generate a model here, test it against the data, compute residuals, ultimately a log-likelihood
	### use a dictionary!
	### if run_LUNA takes (tau0, period, Rplan), for example, you can feed in a dictionary of values
	### with the keys equal to the keywords and it will interpret it! use run_LUNA(**dictionary)

	param_dict = {} ### initialize it
	#for pidx, parlab in enumerate(param_labels):
	for pidx, parlab in enumerate(mn_param_labels):
		param_dict[parlab] = cube[pidx]

	### now you should be able to run_LUNA(param_dict)
	LUNA_times, LUNA_fluxes = pyluna.run_LUNA(data_times, **param_dict, add_noise='n', show_plots='n')

	loglikelihood = np.nansum(-0.5 * ((LUNA_fluxes - data_fluxes) / data_errors)**2) ### SHOULD MAKE THIS BETTER, to super-penalize running out of bounds!
	return loglikelihood 



def mp_multinest(times, fluxes, errors, param_labels, param_prior_forms, param_limit_tuple, nlive, targetID, LUNAfit='y', show_plot='y'):
	### this function will start PyMultiNest!
	"""
	- param_labels is just an array of labels of the parameters you're fitting. MUST MATCH pyluna keywords! 
	- param_prior_forms is an array of prior formats, either 'uniform', 'loguniform', 'normal', 'beta', or 'truncnorm'.
	- param_limit_tuble is an array of tuples dictating the limits or shape parameters of the priors.
	"""

	global mn_param_labels
	global mn_prior_forms
	global mn_limit_tuple
	global data_times
	global data_fluxes
	global data_errors

	mn_param_labels = []
	mn_prior_forms = []
	mn_limit_tuple = []
	data_times = []
	data_fluxes = []
	data_errors = []


	for parlab, parprior, partup in zip(param_labels, param_prior_forms, param_limit_tuple):
		mn_param_labels.append(parlab)
		mn_prior_forms.append(parprior)
		mn_limit_tuple.append(partup)

	for t,f,e in zip(times, fluxes, errors):
		data_times.append(t)
		data_fluxes.append(f)
		data_errors.append(e)


	if os.path.exists('MultiNest_fits'):
		pass
	else:
		os.system('mkdir MultiNest_fits')


	if LUNAfit == 'y':
		outputdir = 'MultiNest_fits/'+str(targetID)
		if os.path.exists(outputdir):
			pass
		else:
			os.system('mkdir '+outputdir)
			os.system('mkdir '+outputdir+'/chains')
		outputdir = outputdir+'/chains'

		pymultinest.run(LogLikelihood=pymn_loglike_LUNA, Prior=pymn_prior, n_dims=len(param_labels), n_live_points=nlive, outputfiles_basename=outputdir+'/'+str(targetID), resume=True, verbose=True)
		json.dump(param_labels, open(outputdir+'/'+str(targetID)+"_params.json", 'w')) ### save parameter names

		"""
		if show_plot=='y':
			try:
				plot.figure()
				plt.plot(times, data, c='r', marker='+')
				a = pymultinest.Analyzer(outputfiles_basename=outputdir+'/'+str(targetID), n_params=len(param_labels))
				
				plot_dict = {}
				for parlab, dict_val in zip(param_labels, a.get_equal_weighted_posterior()[::100,:-1]):
					plot_dict[parlab] = dict_val

				plt.plot(data_times, run_LUNA(data_times, **plot_dict, add_noise='n'), c='b', alpha=0.3)
				plt.show()
			except:
				print("could not plot the solution.")
		"""

def mp_emcee():
	print('nothing happening right now.')

