from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

plt.rcParams["font.family"] = 'serif'




##### this code will put in one place all your plotting choices


def mp_scatter(xarr, yarr, xerr=None, yerr=None, xlabel=None, ylabel=None, label=None, title=None, xscale='linear', yscale='linear', facecolor='LightCoral', edgecolor='k', size=None, colorarr=None, colormap=None, alpha=1, plot_legend='y', show_plot='y'):
	
	#### BASIC SCATTER PLOT
	plt.figure(figsize=(6,8))
	### define your colors
	if type(colormap) != type(None):
		scatter_cmap = cm.get_cmap(colormap)
	else:
		scatter_cmap = cm.get_cmap('viridis')

	if type(colorarr) != type(None):
		norm_colorarr = (colorarr - np.nanmin(colorarr)) / (np.nanmax(colorarr) - np.nanmin(colorarr))
		colors = scatter_cmap(norm_colorarr)
		plt.scatter(xarr, yarr, color=colors, s=size, edgecolor=edgecolor, alpha=alpha, zorder=1, label=None)
	else:
		plt.scatter(xarr, yarr, color=facecolor, edgecolor=edgecolor, s=size, alpha=alpha, zorder=1, label=None)

	plt.errorbar(xarr, yarr, xerr=xerr, ecolor='k', alpha=0.5, fmt='none', zorder=0)
	plt.errorbar(xarr, yarr, yerr=yerr, ecolor='k', alpha=0.5, fmt='none', zorder=0)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.xscale(xscale)
	plt.yscale(yscale)
	if plot_legend == 'y':
		plt.legend()
	plt.title(title)

	if show_plot == 'y':
		plt.show()


def mp_plot(xarr, yarr, xlabel=None, ylabel=None, xscale='linear', yscale='linear', label=None, title=None, color='DodgerBlue', linewidth=1, linestyle='solid', alpha=1, plot_legend='y', show_plot='y'):
	plt.plot(xarr, yarr, color=color, linestyle=linestyle, linewidth=linewidth)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.xscale(xscale)
	plt.yscale(yscale)
	plt.title(title)
	if plot_legend == 'y':
		plt.legend()
	if show_plot == 'y':
		plt.show()




def mp_hist(vals, nbins=20, bins=None, xscale='linear', title=None, yscale='linear', facecolor='DodgerBlue', edgecolor='k', alpha=0.7, xlabel=None, ylabel=None):
	if type(bins) == type(None):
		if xscale == 'linear':
			bins =np.linspace(np.nanmin(vals), np.nanmax(vals), nbins)
		elif xscale == 'log':
			bins = np.logspace(np.log10(np.nanmin(vals)), np.log10(np.nanmax(vals)), nbins)
	plt.figure(figsize=(6,8))
	n,b,e = plt.hist(vals, bins=bins, facecolor=facecolor, edgecolor=edgecolor, alpha=alpha)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.xscale(xscale)
	plt.yscale(yscale)
	plt.title(title)
	plt.show()
	return n,b,e



def mp_heatmap(xarr, yarr, nbins=[20,20], xbins=None, ybins=None, xscale='linear', title=None, yscale='linear', facecolor='DodgerBlue', edgecolor='k', alpha=0.7, xlabel=None, ylabel=None, colormap='viridis'):
	if type(xbins) == type(None) and type(ybins) == type(None):
		if xscale == 'linear':
			xbins =np.linspace(np.nanmin(xarr), np.nanmax(xarr), nbins[0])
			ybins = np.linspace(np.nanmin(yarr), np.nanmax(yarr), nbins[1])
		elif xscale == 'log':
			xbins = np.logspace(np.log10(np.nanmin(xarr)), np.log10(np.nanmax(xarr)), nbins[0])
			ybins = np.logspace(np.log10(np.nanmin(yarr)), np.log10(np.nanmax(yarr)), nbins[1])

	plt.figure(figsize=(8,8))
	h,xe,ye,im = plt.hist2d(xarr, yarr, bins=[xbins,ybins], cmap=colormap)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.xscale(xscale)
	plt.yscale(yscale)
	plt.title(title)
	plt.show()	
	return h,xe,ye,im 









