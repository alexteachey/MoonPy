from __future__ import division
import numpy as np
import astropy
from scipy.stats import norm,beta,truncnorm
import json
import matplotlib.pyplot as plt 
import os
import sys
from mp_batman import run_batman
import pyluna


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

	for pidx, parlabs, parprior, partuple in zip(np.arange(0,len(mn_prior_forms),1), mn_param_labels, mn_prior_forms, mn_limit_tuple):
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



def emcee_lnprior(params):
	#print('NEW CALL.')
	#LUNA_times, LUNA_fluxes = pyluna.run_LUNA(data_times, **param_dict, add_noise='n', show_plots='n')
	in_bounds = 'y'
	#print("parameters = ", params)
	#print(' ')
	#print('mn_limit_tuple = ', mn_limit_tuple)
	#print(' ')
	for param, parlab, parlim in zip(params, mn_param_labels, mn_limit_tuple):
		#print("param name = ", parlab)
		#print('param = ', param)
		#print ('parlim ', parlim)
		#print(" ")
		### test that they are in bounds
		if (param < parlim[0]) or (param > parlim[1]):
			#print("out of bounds!")
			in_bounds = 'n'
			break
	if in_bounds == 'n':
		return -np.inf 
	else:
		return 0.0 


def emcee_lnlike_LUNA(params, data_times, data_fluxes, data_errors):
	param_dict = {} ### initialize it

	for pidx, parlab in enumerate(mn_param_labels):
		param_dict[parlab] = params[pidx]

	LUNA_times, LUNA_fluxes = pyluna.run_LUNA(data_times, **param_dict, add_noise='n', show_plots='n')

	loglikelihood = np.nansum(-0.5 * ((LUNA_fluxes - data_fluxes) / data_errors)**2) ### SHOULD MAKE THIS BETTER, to super-penalize running out of bounds!
	#if (np.isfinite(loglikelihood) == False) or type(loglikelihood) != float:
		#raise Exception('emcee_lnlike_LUNA loglikelihood value is invalid."')

	return loglikelihood 


def emcee_lnlike_batman(params, data_times, data_fluxes, data_errors):
	param_dict = {} ### initialize it

	for pidx, parlab in enumerate(mn_param_labels):
		param_dict[parlab] = params[pidx]

	batman_times, batman_fluxes = run_batman(data_times, **param_dict, add_noise='n', show_plots='n')

	loglikelihood = np.nansum(-0.5 * ((batman_fluxes - data_fluxes) / data_errors)**2) ### SHOULD MAKE THIS BETTER, to super-penalize running out of bounds!
	#if (np.isfinite(loglikelihood) == False) or type(loglikelihood) != float:
	#	print("loglikelihood = ", loglikelihood)
	#	raise Exception('emcee_lnlike_batman loglikelihood value is invalid."')
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




def pymn_loglike_LUNA(cube, ndim, nparams):
	### have to generate a model here, test it against the data, compute residuals, ultimately a log-likelihood
	### use a dictionary!
	### if run_LUNA takes (tau0, period, Rplan), for example, you can feed in a dictionary of values
	### with the keys equal to the keywords and it will interpret it! use run_LUNA(**dictionary)

	param_dict = {} ### initialize it

	for pidx, parlab in enumerate(mn_param_labels):
		param_dict[parlab] = cube[pidx]

	### now you should be able to run_LUNA(param_dict)
	LUNA_times, LUNA_fluxes = pyluna.run_LUNA(data_times, **param_dict, add_noise='n', show_plots='n')

	loglikelihood = np.nansum(-0.5 * ((LUNA_fluxes - data_fluxes) / data_errors)**2) ### SHOULD MAKE THIS BETTER, to super-penalize running out of bounds!
	return loglikelihood 



def pymn_loglike_batman(cube, ndim, nparams):
	### have to generate a model here, test it against the data, compute residuals, ultimately a log-likelihood
	### use a dictionary!
	### if run_LUNA takes (tau0, period, Rplan), for example, you can feed in a dictionary of values
	### with the keys equal to the keywords and it will interpret it! use run_LUNA(**dictionary)

	param_dict = {} ### initialize it

	for pidx, parlab in enumerate(mn_param_labels):
		param_dict[parlab] = cube[pidx]

	### now you should be able to run_LUNA(param_dict)
	batman_times, batman_fluxes = run_batman(data_times, **param_dict, add_noise='n', show_plots='n')

	loglikelihood = np.nansum(-0.5 * ((batman_fluxes - data_fluxes) / data_errors)**2) ### SHOULD MAKE THIS BETTER, to super-penalize running out of bounds!
	return loglikelihood 



def mp_multinest(times, fluxes, errors, param_labels, param_prior_forms, param_limit_tuple, nlive, targetID, modelcode="LUNA", show_plot='y'):
	import pymultinest
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


	if modelcode=='LUNA':
		outputdir = 'MultiNest_fits/LUNA/'+str(targetID)
	elif modelcode=='batman':
		outputdir = 'MultiNest_fits/batman/'+str(targetID)

	if os.path.exists(outputdir):
		pass
	else:
		os.system('mkdir '+outputdir)
		os.system('mkdir '+outputdir+'/chains')
	outputdir = outputdir+'/chains'

	if modelcode == 'LUNA':
		n_params = len(mn_param_labels)
		pymultinest.run(LogLikelihood=pymn_loglike_LUNA, Prior=pymn_prior, n_dims=n_params, n_live_points=nlive, outputfiles_basename=outputdir+'/'+str(targetID), resume=True, verbose=True)
	elif modelcode == "batman":
		n_params = len(mn_param_labels)
		pymultinest.run(LogLikelihood=pymn_loglike_batman, Prior=pymn_prior, n_dims=n_params, n_live_points=nlive, outputfiles_basename=outputdir+'/'+str(targetID), resume=True, verbose=True)
	
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

def mp_emcee(times, fluxes, errors, param_labels, param_prior_forms, param_limit_tuple, nwalkers, nsteps, targetID, resume=True, modelcode="LUNA", storechain=False, burnin_pct=0.1, show_plot='y'):
	import emcee
	import corner

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

	if os.path.exists('emcee_fits'):
		pass
	else:
		os.system('mkdir emcee_fits')

	if os.path.exists('emcee_fits/LUNA'):
		pass
	else:
		os.system('mkdir emcee_fits/LUNA')

	if os.path.exists('emcee_fits/batman'):
		pass
	else:
		os.system('mkdir emcee_fits/batman')

	if modelcode=='LUNA':
		outputdir = 'emcee_fits/LUNA/'+str(targetID)
	elif modelcode=='batman':
		outputdir = 'emcee_fits/batman/'+str(targetID)

	if os.path.exists(outputdir):
		pass
	else:
		os.system('mkdir '+outputdir)
		os.system('mkdir '+outputdir+'/chains')
	outputdir = outputdir+'/chains'
	chainfile = outputdir+'/'+str(targetID)+'_mcmc_chain.txt'
	lnpostfile = outputdir+'/'+str(targetID)+"_mcmc_lnpost.txt"


	### if chainfile doesn't exist, resume == False!
	if os.path.exists(chainfile) == False:
		resume = False 

	### INITIALIZE EMCEE PARAMETERS
	ndim = len(mn_param_labels)
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
			for partup in mn_limit_tuple:
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
		sampler = emcee.EnsembleSampler(nwalkers, ndim, emcee_lnprob_LUNA, args=(data_times, data_fluxes, data_errors))

	elif modelcode == "batman":
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
	#    if (i+1) % 100 == 0:
	#       print("{0:5.1%}".format(float(i) / nsteps))
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
		fig = corner.corner(samples, labels=mn_param_labels)
		plt.savefig(outputdir+'/'+str(targetID)+"_corner.png")
		plt.close()


