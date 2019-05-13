### Welcome to MoonPy!
README last updated: May 13, 2019.

This document will walk you through the basics of using the MoonPy code. 
In time this code will become more sophisticated, but right now we can just do
a few things.

1.) INITIALIZE A LIGHT CURVE OBJECT. Proper usage is:

>>> lc_objectname = MoonPy_LC(lc_times=None, lc_fluxes=None, lc_errors=None, targetID=None, target_type='koi', quarters='all', telescope='kepler', RA=None, Dec=None, coord_format='degrees', search_radius=5, lc_format='pdc', sc=False, ffi='y', lc_meta=None, save_lc='y', loadfile='n')

This object is designed to be versatile. You can either
a) supply times, fluxes, and errors as arrays;
b) supply a targetID (either a KOI, Kepler planet, or KIC); or
c) supply coordinates for an object search

If you choose option (b), you must specify the target_type as "koi", "planet", or "kic". 
If it's a 'koi', input the targetID as a the koi number including the decimal as a float. Example: 5084.01
If it's a 'planet', input the targetID as the planet number and 'b' as a string. Example: "1625b".
If it's a 'kic', input the targetID as the KIC number without the prefix, as a float. Example: 4760478

The coordinate search (c) performs a cone search with a 5 arcsecond radius through Simbad. You may change the 
cone size by adjusting the "search_radius" keyword. Some targets have multiple aliases, and if the first hit
is not either a KOI, kepler planet or KIC star, an attempt will be made to find this name amongst the aliases.
Also note that your options for coord_format are 'degrees' and "sexagesimal", but if you input sexagesimal
without indicating it an attempt is made to recognize this and change the coord_format on the fly. Appropriate
syntax for sexagesimal coordinates is '9h36m43.5s' for RA and '+39d42m46.83s' for Dec. I think you can go arbitrarily
precise with the decimal places but I haven't tested this extensively.

You may also download only select quarters if you wish by supplying an array of quarter numbers in the 'quarters' keyword.
Currently supported formats for the light curve download is "sap" and "pdc". If you wish to download both
(possibly included in a future release) you should just initialize two different light curve objects.

The first time you use the code you will want to specify your SAVEPATH at the very top of the code.
Unless you specify otherwise your light curve will be saved in the savepath ('save_lc' keyword).

Other functionality listed above is forthcoming.

##########


2.) PLOT THE DATA.

Plotting the data is simple. Once you have generated your light curve object (step 1 above), you can
plot the light curve simply by calling 

>>> lc_objectname.plot(facecolor='LightCoral', edgecolor='k', errorbar='n', quarters='all', include_flagged='n', detrended='y'))

If the light curve has already been detrended, you will see the detrended light curve. IF NOT, you will get a 
warning that the light curve has not yet been detrended and you will see instead the raw light curve.


3.) DETREND THE DATA.

At this time the only detrending option is the Cosine Filtering Autocorrelation Minimization (CoFiAM) algorithm.
This is the default. Future options may include Dan Foreman-Mackey's "untrendy" and "george" packages, and possibly
others. The usage is simple:

>>> lc_objectname.detrend(dmeth='cofiam', save_lc='y')

As before, you may choose to save this light curve or not by altering the "save_lc" keyword. Once again,
you will need to fill in your savepath the first time you use the code.



4.) FUTURE FUNCTIONALITY

Additional detrending options are in the works, as are implementations of MultiNest and emcee. BATMAN support
is also planned, and (hopefully), implementation of David Kipping's LUNA code for generating planet+moon light curves.
STAY TUNED!










