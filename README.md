![MoonPy_logo](https://github.com/alexteachey/MoonPy/blob/master/moonpy_logo.png)



# Welcome to MoonPy!
last update: 9 August 2022. Developed by [Alex Teachey](http://www.alexteachey.com).

This document will walk you through the basics of using the MoonPy code. 

**MoonPy is designed to make downloading, plotting, detrending, and *fitting PLANET and MOON models to light curves a breeze.*** These are routine tasks, but often require a lot of coding to get going with. Sometimes you just want to take a quick look and have all the numbers at your fingertips -- enter MoonPy!

Of course, the easier a tool is to use, the more control you are handing over to the code / programmer. As such, *MoonPy is intended mainly as a first analysis tool, not necessarily a final analysis tool.* Please keep this in mind!

MoonPy simply would not exist if not for a number of incredibly powerful packages that have been developed by other scientists (more on that below). Thank you to all those developers! And **if you use MoonPy in your research please remember to cite/acknowledge the authors of the packages that have been utilized here. We'd appreciate a shoutout, too... please cite:**

[\bibitem[Teachey \& Kipping(2021)]{2021MNRAS.508.2620T} Teachey, A. \& Kipping, D.\ 2021, \mnras, 508, 2620. doi:10.1093/mnras/stab2694](https://ui.adsabs.harvard.edu/abs/2021MNRAS.508.2620T/abstract).

 [![DOI](https://zenodo.org/badge/186446465.svg)](https://zenodo.org/badge/latestdoi/186446465)

A special shoutout is due to Hippke & Heller for introducing [the ```Pandora``` code](https://github.com/hippke/Pandora), as well as Gordon & Agol for introducing [the ```gefera``` code.](https://github.com/tagordon) These are the first publicly-available codes for modeling exomoon light curves, which is a big deal. The original aim of the ```MoonPy``` package was to make exomoon modeling widely accessible, and it is now possible to realize that objective thanks to their codes. MoonPy is providing a wrapper that makes it even easier to implement their already user-friendly codes. 

 

Let's walk through the basics.



# Installation
Installation is now easier than ever. Begin by cloning the GitHub repository wherever you see fit. For example, within the desired directory, call:

```>>> git clone https://github.com/alexteachey/MoonPy```

Then navigate to the MoonPy directory you have just cloned, and within the terminal type

```python install_moonpy.py```

This will 1) create a new conda environment for MoonPy with a name of your choosing, 2) install a variety of required and recommended packages within that environment, and 3) create the path so that your Python can see MoonPy as an Python-importable package. After this is complete, simply 

```conda activate ENV_NAME```

where ```ENV_NAME``` is the name of the conda environment you have just established. You should now be able to type

```from moonpy import *```

And you're off to the races. **Please note:** The first time you boot up MoonPy it may take a few moments to load. After that it should be much faster!

*If the installation process above fails, dependencies below can be installed individually.*

* **LINUX USERS:** To ensure MoonPy loads correctly, users may need to call* ```python pathmaker.py``` *within the MoonPy directory, and with the desired conda environment activated. If you opt to install VESPA (recommended), you may also need to call* ```python vespa_script_updater.py``` *to fix a deprecated keyword in some vespa scripts. This should be done within the separate VESPA conda environment. (breakdowns include a deprecated 'representation' keyword, and/or a missing 'TESS' keyword. These errors are fixed by vespa_script_updater.py .*

**A list of package dependencies is listed at the bottom of this page**.


## NB: Frequent Updates!
Be advised that MoonPy is being updated all the time, and if I'm being honest, I'm not always adhering to best practices in terms of maintaing stable and beta versions. By using MoonPy you acknowledge that you are on this ride with me... things may break from time to time, and if they do, I would appreciate hearing about them. This also means, it's important to keep MoonPy up-to date by cloning the latest version of the repository.



# OVERVIEW
Below you will find the details of some of these tools. At present the user will want to utilize the following functionality:

* **Initialize a light curve object** using ```MoonpyLC()```
* **Access adopted NASA Exoplanet Archive attributes** using [the same column names from the archive](https://exoplanetarchive.ipac.caltech.edu/docs/API_PS_columns.html).
* **plot the light curve** using ```plot_lc()```.
* **generate a Lomb-Scargle** periodogram using ```genLS()```.
* **Detrend the light curve** using ```detrend()```
* **fit a planet transit model** (using the ```exoplanet```package) to the light curve using ```run_planet_fit()```.
* **fit PLANET and MOON models** (using ```Pandora``` or ```gefera```, and ```UltraNest```) using ```fit(modelcode='Pandora', model='M')```.
* **compute the Bayes factor** to find evidence for/against the moon using ```moon_evidence()```.
* **run a false positive calculation with VESPA** using ```run_vespa()```. ***Must activate separate vespa environment.***
* **access a variety of handy tools** in the ```mp_tools.py``` module.

*And more!*



# QUICKSTART

## Initialize a light curve object.
MoonPy is built around manipulating light curves as objects. To start, simply import moonpy:

```>>> from moonpy import *```

Then you can initialize an object:

```>>> lc_object=MoonpyLC(targetID=TARGET_NAME, clobber='n')```

targetID can take a range of common names for stars / planets(e.g. 'Kepler-1625b', 'KOI-5084.01', 'KIC4760478', 'TOI-216.01', 'TIC55652896', 'EPIC201170410'...)

MoonPy will try to infer which telescope's data should be accessed based on the name of the object (KICs, KOIs, and Kepler- planets get Kepler data, K2 or EPIC numbers retrieve K2 data, TOIs and TICs access TESS data), but you may also specify the telescope if the target name is not one of these, or if you want to search for data from a different telescope. For example, some Kepler planets were also observed by TESS, so simply specify *telescope='tess'*, etc, for more controlled access.

It is also possible to use MoonPy tools on a light curve that the user supplies. In this case, the user will want to supply arrays for ```lc_times```, ```lc_fluxes```, and ```lc_errors``` when initializing a MoonpyLC object. In addition, the ```usr_dict``` containing important parameters for the planet can be included (keys are "period", "tau0", "impact", "duration_hours", "rprstar", "sma_AU", and "rp_rearth"). If these are not supplied, the user will be prompted to enter them.




## Visualizing your data
You can then plot the light curve with

```>>> lc_object.plot_lc()```

A plot of all the available data should appear, and known transits should be flagged. If the light curve has already been detrended (see below), you will see two plots, the top one being the raw data with the trend line overplotted, the bottom one showing the detrended light curve. If planet and / or moon models have been produced, they are also overplotted by default.

You can also generate a Lomb-Scargle periodogram simply by calling 

```>>> lc_object.genLS()```

This produces quarter-by-quarter plots that may indicate the presence of periodic signals (likely rotational information about the star).

After you have run a planet or moon model (see below), you can generate an animation of that model (using median values from the posteriors) using

```>>> lc_object.animate_moon()```.




## Data Manipulation

To detrend a light curve, simply call

```>>> lc_object.detrend(dmeth='cofiam')```

A variety of detrending methods are available ('cofiam', 'median_value', 'phasma', 'polyAM', 'polyLOC', moving_median', and 'methmarg'). **Users are strongly advised to examine the results of the detrending before using for science!** 

Then you can ```plot_lc()``` again to see both the original data, the trend model, and the detrended light curve.




## Accessing Planet Attributes 

You can see the full list of attributes associated with this object (values, arrays, and functions) by calling

```>>> lc_object.print_attributes()```

Many of these [attributes are taken directly from the NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/docs/API_PS_columns.html), and have the same name. Therefore you can call up many values of interest for this target simply by calling the attribute of choice, for example:

```>>> k1625.pl_orbper```

```287.378949```

MoonPy also automatically queries the [SIMBAD astronomical database](https://simbad.unistra.fr/simbad/) to identify aliases for the planet/star in question. You can see these with 

```>>> lc_object.find_aliases()```.

MoonPy also automatically identifies all other confirmed planets in the system, and *should* flag them when plotted. To see any neighbors, you can call

```>>> k1625.find_neighbors()```

You may call ```lc_object.get_neighbors()```, which will generate a dictionary of neighbor light curve objects, each with their own attributes. For example

```>>> k167e.neighbor_dict```

returns

```{'K167b': <moonpy.MoonpyLC object at 0x7f80a5ff17f0>, 'K167c': <moonpy.MoonpyLC object at 0x7f80a3ab87c0>, 'K167d': <moonpy.MoonpyLC object at 0x7f80a3ab89a0>}```




## Model Fitting 
From here, you might opt to run a planet transit fit with ```Exoplanet```. Simply call

```>>> lc_object.run_planet_fit()```

And you're off to the races. These fits tend to take a few minutes, and jointly fit a planet transit model with a GP trend model.

For a more in-depth analysis, you can **fit planet and moon models using ```Pandora``` and ```UltraNest```.** These fits may take a few hours. ```MoonPy``` retrieves established NASA Exoplanet Archive attributes for these planets and uses them as Gaussian priors for your model fitting. Simply call

```>>> lc_object.fit(modelcode='Pandora', model='M')``` 

to start your run. You can also set ```modelcode='gefera'``` and set ```model='P'``` to run a 'P'lanet model instead of a 'M'oon model.

After both models have run, you may call 

```>>> lc_object.moon_evidence()```

Which will access the Bayesian evidence for each model, compare them, and let you know whether the planet or moon model is favored, and to what degree.

You can also **access the model posteriors** by calling 

```>>> lc_object.get_Pandora_posteriors(model='M') # or P```

which will load the model posteriors as attributes for the ```lc_object```. Then call 

```>>> lc_object.Pandora_moon_PEWdict```

or

```>>> lc_object.Pandora_planet_PEWdict```

to retrieve dictionaries for the two models. (NB: PEW standards for "Posterior Equal Weights").




## False Positive Vetting

There are a variety of astrophysical scenarios that can mimic a transiting planet. In particular, eclipsing binaries are a major source of false positives. Tim Morton's ```VESPA``` code was designed to compute a number of false positive probabilities.

Unfortunately, ```VESPA``` is now somewhat tricky to get installed and running, but MoonPy has you covered!. **Ensure that during installation you opted to install VESPA, which places VESPA in its own conda environment.** Due to older dependencies it is important that ```VESPA``` be quarantined in this way. Assuming this has gone through successfully, you can activate your ```VESPA``` specific conda environment, and then call 


```>>> lc_object.run_vespa()```

```MoonPy``` automatically generates the you will generate the necessary input files for this target. 




## Data storage

Each target gets its own directory within the ```Central_Data``` directory, which by default is saved within the MoonPy directory. These directories are separated as Kepler, K2, and TESS. You can see the location of these files by calling ```lc_object.savepath```. VESPA runs will also be carried out and saved within this same directory, and it *should* be possible to run multiple VESPA fits at once as a result.




## Up-To-Date Data
New TESS data continues to arrive, so **MoonPy automatically identifies when new sectors come online and will download them if available**. The planet archives (NASA Exoplanet Archive, and ExoFOP) are also regularly updated, so MoonPy checks to see how old your databases are and gives the option (once per day) of downloading a new version. Be advised, this can sometimes take several minutes, so if you are looking at a well-known target, chances are the parameters are more-or-less fixed.




## Handy Tools
MoonPy utilizes a variety of functions behind the scenes that exoplanet astronomers frequently need. So they might come in handy. Here are some of them (please sanity check the results before using them!):

* ```effective_radius(density, mass)```
* ```RHill(sma_plan, m_star, m_plan)```
* ```Kep3_pfroma(sma, m1, m2, sma_unit = 'meters', output_format='days', val_only='y')```
* ```Kep3_afromp(period, m1, m2, val_only='y', unit='days')```
* ```mass_from_density(density, radius)```
* ```inc_from_impact(impact, rstar, sma, unit='radians')```
* ```impact_from_inc(inclination, rstar, sma, unit='degrees')```
* ```Tdur(period, Rstar, Rplan, impact, sma)``` # assumes circular orbit
* ```deg2rad(degrees)```
* ```rad2deg(radians)```
* ```Roche(Rsat, Mplan, Msat)```
* ```density_from_orbit(a_over_R, Porbit, in_unit='days', out_unit='mks')```
* ```lc_fold(times, fluxes, errors, tau0, period, phase_offset=0.0)```
* ```DWstat(data, model)```
* ```flux_from_mags(target_mag, ref_mag, ref_flux)```




## Dependencies

This code **requires** the following standard packages, which *should* be installed following the instructions above! Dependencies get *very* tricky very fast, so by far the best option is to use the MoonPy installer, or if you have issues there, create a conda environment from one of the included .yml files.

If users have installation issues, they may try to install individually using the instructions below:

* [astropy](https://www.astropy.org/) -- conda install -c anaconda astropy
* [astroquery](https://astroquery.readthedocs.io/en/latest/) -- conda install c astropy astroquery
* [matplotlib](https://matplotlib.org/) -- conda install -c conda-forge matplotlib
* [pandas](https://pandas.pydata.org/) -- conda install -c anaconda pandas
* [scipy](https://scipy.org/) -- conda install -c anaconda scipy
* [numba](https://numba.pydata.org/) -- conda install -c conda-forge numba

For full functionality, the following packages will also be needed (but MoonPy should be able to be run without them):
* [exoplanet](https://docs.exoplanet.codes/en/latest/) -- conda install -c conda-forge exoplanet
* [pandora](https://github.com/hippke/Pandora) -- pip install pandoramoon
* [gefera](https://github.com/tagordon/gefera) -- pip install gefera 
* [pymc3](https://docs.pymc.io/en/v3/) -- conda install -c conda-forge pymc3
* ```pymc3_ext``` -- conda install -c conda-forge pymc3_ext
* [arviz](https://arviz-devs.github.io/arviz/) -- conda install -c conda-forge arviz=0.11.0 **NOTE:** this version is important to avoid a breakdown with ```pymc3```.
* [corner](https://github.com/dfm/corner.py) -- conda install -c astropy corner
* [celerite2](https://celerite2.readthedocs.io/en/latest/) -- conda install -c conda-forge celerite2
* ```aesara_theano_fallback``` -- conda install -c conda-forge aesara-theano-fallback
* [george](https://george.readthedocs.io/en/latest/) -- conda install -c conda-forge george 
* [untrendy](https://github.com/dfm/untrendy) -- pip install untrendy
* [emcee](http://dfm.io/emcee/current/) -- conda install -c conda-forge emcee
* [batman](https://www.cfa.harvard.edu/~lkreidberg/batman/) - pip install batman-package
* [PyMultiNest](https://johannesbuchner.github.io/PyMultiNest/) - pip install pymultinest 
* [UltraNest](https://johannesbuchner.github.io/UltraNest/index.html) -- pip install ultranest
* [vespa](https://github.com/timothydmorton/VESPA) -- recommended to install via MoonPy installer
* [isochrones](https://isochrones.readthedocs.io/en/latest/) -- recommended to install via MoonPy installer


Package imports that are not standard and are not included in the moonpy package (i.e. the packages above) should all be imported within the relevant function that utilizes them. Therefore, you can hopefully boot up moonpy and use it even if you lack some of the above distributions.
