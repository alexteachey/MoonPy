from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mp_tools import * 

plt.rcParams["font.family"] = 'serif'




##### this code will put in one place all your plotting choices





def mp_scatter(xarr, yarr, xerr=None, yerr=None, xlabel=None, ylabel=None, label=None, title=None, xscale='linear', yscale='linear', facecolor='LightCoral', edgecolor='k', size=None, colorarr=None, colormap=None, alpha=1, zorder=1, plot_legend='y', show_plot='y'):
	
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
		plt.scatter(xarr, yarr, color=colors, s=size, edgecolor=edgecolor, alpha=alpha, zorder=zorder, label=None)
	else:
		plt.scatter(xarr, yarr, color=facecolor, edgecolor=edgecolor, s=size, alpha=alpha, zorder=zorder, label=None)

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


def mp_scatter_and_hist(x, y, xlabel=None, ylabel=None, xline=None, yline=None, xscale='log', yscale='log', facecolor='LightCoral'):
	# definitions for the axes
	left, width = 0.15, 0.60
	bottom, height = 0.1, 0.65
	spacing = 0.03

	rect_scatter = [left, bottom, width, height]
	rect_histx = [left, bottom + height + spacing, width, 0.2]
	rect_histy = [left + width + spacing, bottom, 0.2, height]

	# start with a rectangular Figure
	plt.figure(figsize=(8, 8))

	ax_scatter = plt.axes(rect_scatter)
	ax_scatter.tick_params(direction='in', top=True, right=True)
	ax_histx = plt.axes(rect_histx)
	ax_histx.tick_params(direction='in', labelbottom=False)
	ax_histy = plt.axes(rect_histy)
	ax_histy.tick_params(direction='in', labelleft=False)

	xmin, xmax = np.nanmin(x), np.nanmax(x)
	ymin, ymax = np.nanmin(y), np.nanmax(y)

	if xscale == 'log':
		ax_scatter.set_xscale('log')
		ax_histx.set_xscale('log') ### may need to be ax_histx
		ax_scatter.set_xlim(xmin, xmax)
		ax_histx.set_xlim(xmin, xmax)

	if yscale == 'log':
		ax_scatter.set_yscale('log')
		ax_histy.set_yscale('log') #### may need to be ax_histy
		ax_scatter.set_ylim(ymin, ymax)
		ax_histy.set_ylim(ymin, ymax)

	# the scatter plot:
	ax_scatter.scatter(x, y, facecolor=facecolor, edgecolor='k', s=20, alpha=0.7)

	# now determine nice limits by hand:
	binwidth = 0.25
	#lim = np.ceil(np.abs([x, y]).max() / binwidth) * binwidth
	#ax_scatter.set_xlim((-lim, lim))
	#ax_scatter.set_ylim((-lim, lim))
	ax_scatter.set_xlim(np.nanmin(x), np.nanmax(x))
	ax_scatter.set_ylim(np.nanmin(y), np.nanmax(y))
	ax_scatter.set_xlabel(xlabel)
	ax_scatter.set_ylabel(ylabel)

	ax_scatter.tick_params(axis='x', labelsize=16)
	ax_scatter.tick_params(axis='y', labelsize=16)

	if xline != None:
		xline_xvals = np.linspace(xline, xline, 100)
		xline_yvals = np.linspace(np.nanmin(y), np.nanmax(y), 100)
		ax_scatter.plot(xline_xvals, xline_yvals, c='k', linestyle='--')
	if yline != None:
		yline_xvals = np.linspace(np.nanmin(x), np.nanmax(x), 100)
		yline_yvals = np.linspace(yline, yline, 100)
		ax_scatter.plot(yline_xvals, yline_yvals, c='k', linestyle='--')

	if xlabel != None:
		ax_scatter.set_xlabel(xlabel, fontsize=20)
	if ylabel != None:
		ax_scatter.set_ylabel(ylabel, fontsize=20)

	#bins = np.arange(-lim, lim + binwidth, binwidth)
	#bins = 20
	#xbins = np.logspace(np.log10(np.nanmin(x)), np.log10(np.nanmax(x)), 20)
	#ybins = np.logspace(np.log10(np.nanmin(y)), np.log10(np.nanmax(y)), 20)
	#xbins = np.logspace(-3,5,30)
	#ybins = np.logspace(-3,5,30)
	xbins = np.logspace(np.log10(xmin), np.log10(xmax), 30)
	ybins = np.logspace(np.log10(ymin), np.log10(ymax), 30)

	ax_histx.hist(x, bins=xbins, facecolor=facecolor, edgecolor='k', alpha=0.7)
	ax_histy.hist(y, bins=ybins, orientation='horizontal', facecolor=facecolor, edgecolor='k', alpha=0.7)

	ax_histx.tick_params(axis='y', labelsize=16)
	ax_histy.tick_params(axis='x', labelsize=16)

	#ax_histx.set_xlim(ax_scatter.get_xlim())
	#ax_histy.set_ylim(ax_scatter.get_ylim())

	plt.show()



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









