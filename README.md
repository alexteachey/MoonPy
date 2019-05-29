# Welcome to MoonPy!
last updated: May 29, 2019.

This document will walk you through the basics of using the MoonPy code. 
In time this code will become more sophisticated, but right now we can just do
a few things.


## Prerequisites

This code requires a number of packages, including: 

* [PyMultiNest](https://johannesbuchner.github.io/PyMultiNest/) - standard parameter estimator and Bayesian evidence calculator, developed by Johannes Buchner as a wrapper to the standard MultiNest code (Feroz and Hobson)
* [emcee](http://dfm.io/emcee/current/) -- alternative to MultiNest, for parameter estimation, developed by Dan Foreman-Mackey (DFM) et al
* [batman](https://www.cfa.harvard.edu/~lkreidberg/batman/) - Standard transit-modeling code by L. Kreidberg
* [kplr](http://dfm.io/kplr/) - for downloading Kepler light curves, developed by DFM
* [k2plr](https://github.com/rodluger/k2plr) -- a modification of the kplr code for downloading K2 data, by Rodrigo Luger
* [corner](https://github.com/dfm/corner.py) -- for visualizing fit parameters, developed by DFM
* [george](https://george.readthedocs.io/en/latest/) -- a Gaussian process regression tool developed by (you guessed it) DFM
* [untrendy](https://github.com/dfm/untrendy) -- another DFM package for fast and easy detrending (currently having difficulties here)

As well as standard packages (numpy, scipy, astropy, matplotlib, pandas) that likely came with your python distribution.

If for some reason you have difficulties with one or more of these packages, my best advice is to go into this source code, deactivate
the imports, and proceed without using these particular packages. 


## USAGE:

Eventually this package will come with a pip install, or something along those lines. For now however, you'll just want to clone this repository and run everything locally. If you keep all these scripts together they ought to talk to each other. 

Within the MoonPy directory, simply use 

```
from moonpy import *
```
and that ought to get you going.


## INITIALIZE A LIGHT CURVE OBJECT.

Proper usage is:
```
>>> lc_object = MoonpyLC(targetID=None, lc_times=None, lc_fluxes=None, lc_errors=None, target_type=None, quarters='all', telescope=None, RA=None, Dec=None, coord_format='degrees', search_radius=5, lc_format='pdc', remove_flagged='y', sc=False, ffi='y', save_lc='y', load_lc='n')
```

### KEYWORDS

*targetID*: may be a Kepler planet, a KIC star, a KOI, or a K2 target. If you leave off the prefix you must specify the telescope. Example: "Kepler-1625b", KOI-5084.01, or KIC 4760478. Variations on these (with and without prefix) are (hopefully!) handled to your satisfaction.

*lc_times*: array of times (not yet supported).

*lc_fluxes*: array of fluxes. (not yet supported).

*lc_errors*: array of errors. (not yet supported).

*target_type*: May be "kic", "koi", or "planet" (for confirmed planets a la Kepler-1625b). The code will attempt to intuit this.

*quarters*: acceptable values are "all" or an array of quarters.

*telescope*: may be "Kepler" / "kepler", "K2" / "k2", and eventually "TESS" / "Tess" / "tess".

*RA*: sexagesimal of the form 12h34m14.5s or decimal form.

*Dec*: sexagesimal of the form +50d38m14.25s or decimal form.

*coord_format*: format of your supplied coordinates. An attempt is made to intuit this.

*search_radius*: in arcseconds, the size of the search cone when supploying coordinates.

*lc_format*: May be "sap" (simple aperture photometry) or "pdc" (pre-search data conditioning). Default is "pdc".

*remove_flagged*: automatically removes data points that have quality flags != 0.

*sc*: Boolean, stands for "short cadence". Not doing anything right now.

*ffi*: stands for "Full-Frame Images." Not doing anything right now.

*save_lc*: option to save your light curve once you've generated it as a .csv file.

*load_lc*: if 'y', an attempt is made to load a light curve file you have already generated through this code.


**NOTES:**
This object is designed to be versatile. You can either
a) supply times, fluxes, and errors as arrays (EVENTUALLY)!;
b) supply a targetID (either a KOI, Kepler planet, KIC, or K2 planet) and the name of the telescope; or 
c) supply coordinates for an object search

If you choose option (b), you may need to make it explicit somehow which telescope you want to use.
For example, you can either enter a targetID like "Kepler-1625b", "KOI-5084.01", or "KIC4760478", OR
you may enter "1625b", "5084.01", or "4760478" for the targetID and specify the telescope as "kepler".
The code will do its best to determine the telescope. It should also accept "Kepler" as well as "kepler",
and "K2" as well as "k2". TESS support is in the works. If you have already downloaded this light curve, 
you may set load_lc='y' to attempt to load a file you have already generated (handy if you've already 
detrended the light curve and don't want to do it again.)

The coordinate search (c) performs a cone search with a (default) 5 arcsecond radius through Simbad. You may change the 
cone size by adjusting the "search_radius" keyword. Some targets have multiple aliases, and if the first hit
is not either a KOI, kepler planet or KIC star, an attempt will be made to find this name amongst the aliases.
Also note that your options for coord_format are 'degrees' and "sexagesimal", but if you input sexagesimal
without indicating it an attempt is made to recognize this and change the coord_format on the fly. Appropriate
syntax for sexagesimal coordinates is '9h36m43.5s' for RA and '+39d42m46.83s' for Dec. I think you can go arbitrarily
precise with the decimal places but I haven't tested this extensively. Spaces between hours/minutes/seconds and degrees/minutes/seconds
should be OK.

You may also download only select quarters if you wish by supplying an array of quarter numbers in the 'quarters' keyword.
Currently supported formats for the light curve download is "sap" and "pdc". If you wish to download both
(possibly included in a future release) you should just initialize two different light curve objects.

The first time you use the code you will want to specify your SAVEPATH at the very top of the code.
Unless you specify otherwise your light curve will be saved in the savepath ('save_lc' keyword).

Other functionality listed above is forthcoming, including support for TESS light curves! 




## PLOT THE DATA.

Plotting the data is simple, and I expect the keywords are all self-explanatory. 

Once you have generated your light curve object (step 1 above), you can plot the light curve simply by calling 
```
>>> lc_objectname.plot(facecolor='LightCoral', edgecolor='k', errorbar='n', quarters='all', include_flagged='n', detrended='y'))
```
If the light curve has already been detrended, you will see the detrended light curve. IF NOT, you will get a 
warning that the light curve has not yet been detrended and you will see instead the raw light curve. 

You may also set folded='y' to see a phase fold of the light curve.





## DETREND THE DATA.

At this time the only detrending option is the Cosine Filtering Autocorrelation Minimization (CoFiAM) algorithm "cofiam", and
a median filter "medfilt". Dan Foreman-Mackey's "untrendy" package is also supported, but right now Scipy is throwing
an error saying that the x-values must be strictly increasing. Not sure what that's about. 

Anyway the usage is simple:
```
>>> lc_object.detrend(dmeth='cofiam', save_lc='y', mask_transits='y', skip_ntqs='n', kernel=None, max_degree=30)
```

### KEYWORDS

*dmeth*: currently supported are "cofiam", "medfilt", and "untrendy" (though the last is breaking a lot).
save_lc: default is on. Note that this will overwrite the light curve file you've generated by default in initializing the object, so that now there are five columns: times, fluxes, errors, and fluxes and errors from the detrending. Once again, you will need to fill in your savepath the first time you use the code.

*mask_transits*: by default, transits are masked by calculating the transit times (assuming linear ephemeris), and one full transit duration
is masked on either side of the transit midtime. That is, the total mask is twice the width of the transit. 

*skip_ntqs*: this option allows you to only detrend the quarters that actually contain a transit of the planet you're interested in.
This can be useful with cofiam, for example, since each quarter can take ~1 minute to detrend. Off by default.

*kernel*: this allows you specify the size of the kernel for median filtering.

*max_degree*: for cofiam, this is maximum order k you'll allow cofiam to explore. Practically speaking anything much above this
because too computationally expensive.




## GET PROPERTIES.

New feature as of May 14, 2019 is the "get_properties()" method. This function queries the NASA Exoplanet Archive
to retrieve your target's impact parameter, transit duration, orbital period, all transit midtimes within the 
dataset baseline, the ratio of radii, and the isolated radii for the planet and the star. Uncertainties
for many of these parameters are also available as a tuple, quoting lower and upper sigmas.

The following attributes are supported;
```
lc_object.period # days
lc_object.period_err # tuple
lc_object.tau0 # BKJD
lc_object.tau0_err # tuple
lc_object.impact 
lc_object.impact_err # tuple
lc_object.duration_hours
lc_object.duration_hours_err # tuple 
lc_object.duration_days
lc_object.duration_days_err # tuple
lc_object.rprstar 
lc_object.rprstar_err # tuple
lc_object.rp_rearth # units of Earth radii (native unit on NASA Exoplanet Archive)
lc_object.rp_rearth_err # tuple
lc_object.rp_rjup # units of Jupiter radii (converted without uncertainties)
lc_object.rstar_rsol # units of Solar radii (converted from Rp/Rstar, without propagating uncertainties)
lc_object.depth 
lc_object.taus # all BKJD transit midtimes in the baseline, assuming linear ephemeris
```
If an attribute doesn't come with an uncertainty tuple, it's probably because this is a non-native value
and I haven't bothered to propagate the uncertainties. Will try to implement this in the future.

More attributes may be added in the future. The first time you run this function it will download an ascii table
from NASA Exoplanet Archive, but should not download the table again until 24 hours have elapsed.






## FIT A TRANSIT MODEL TO THE LIGHT CURVE.
(NOTE: the LUNA code developed by D. Kipping is not currently available on GitHub. The pyluna.py script is simply
a wrapper for this code).

You may fit a LUNA or BATMAN model to your detrended data using the following command:

```
>>> lc_object.fit(custom_param_dict=None, fitter='multinest', modelcode='LUNA', skip_ntqs='y', model='M', nlive=500)
```

The *fitter* may be either "multinest" or "emcee" (the latter is still a bit buggy!)

Two *modelcode* options ("LUNA" and "batman") may be used. Right now only one model type "M" is supported for LUNA. 

Future functionality will include evidence testing for four models: P, T, Z, and M.

As used in Teachey & Kipping (2018), the four models are as follows (once they're supported):
(P): a planet-only model that assumes strict linear ephemeris;
(T): a planet-only model that allows the transit times to be fit individually (maximum 6 transits);
(Z): a moon model that sets the moon radius to zero (useful for testing the dynamical effects of the moon); and
(M): a fully physical moon model.

custom_param_dict allows the user to supply their own parameter dictionary to override the default parameter 
values that are supplied in the fitting. The dictionary is as follows:

param_dict['keyword'] = ['prior_type', (prior_lower_bound, prior_upper_bound)]

they keyword must be one of the keywords accepted by pyluna: 
[RpRstar, rhostar, bplan, Pplan, tau0, q1, q2, rhoplan, sat_sma, sat_phase, sat_inc, sat_omega, MsatMp, RsatRp]

the 'prior_type' may be 'uniform', 'loguniform', 'normal', 'lognormal', or 'beta'.

If you wish to keep the default parameters there is no need to supply these.



## PLOT YOUR BEST MODEL..

Once you've run a model fit (either with PyMultiNest or emcee, and with either LUNA or batman) you may wish to overplot your best fit onto the data. This is still under development, but you can use

```
lc_object.plot_bestmodel(self, fitter, modelcode, burnin_pct=0.1)
```
to generate a plot of the resulting model. Note that the "best fit" for each parameter is simply the median value from the chains.
If you want to do something more sophistcated you'll probably want to go in and manually alter the code. 

As of May 29, 2019 this is only working for runs made with emcee. Pymultinest support is coming soon.


## MAKE A CORNER PLOT.

You can also generate a corner plot based on your emcee chains (multinest support coming soon.) Usage here is:

```
lc_object.plot_fit(self, fitter='multinest', modelcode='batman', burnin_pct=0.1)
```

This will save a corner plot in the chains directory where your planet chains were saved.







