import sys 
import os


#### FIX VESPA -- OUT OF DATE KEYWORDS FROM ASTROPY
### '/opt/miniconda3/envs/vespa_with_moonpy_for_mac/lib/python3.6/site-packages/vespa/orbits'
pythonpaths = sys.path 

for pythonpath in pythonpaths:
	if 'site-packages' in pythonpath:
		desired_path = pythonpath
orbitsdir = desired_path+'/vespa/orbits'
starsdir = desired_path+'/vespa/stars'


#### UPDATE POPULATIONS.PY with ORBITS directory
popsfile = open(orbitsdir+'/populations.py', mode='r')
new_popsfile = open(orbitsdir+'/populations_new.py', mode='w')

for nline, line in enumerate(popsfile):
	newline = line.replace('representation=', 'representation_type=')
	new_popsfile.write(newline)

popsfile.close()
new_popsfile.close()
os.system('mv '+orbitsdir+'/populations.py '+orbitsdir+'/populations_DEPRECATED.py')
os.system('mv '+orbitsdir+'/populations_new.py '+orbitsdir+'/populations.py')

print('UPDATED '+orbitsdir+'/populations.py')



#### UPDATE POPULATIONS.PY within STARS directiory
##### THIS IS FOR isochrones 1.1.1

isochrones_version = os.popen('pip freeze | grep isochrones').read()[:-1] ### gets rid of fina \n in the output
print('isochrones_version = ', isochrones_version)
if isochrones_version == 'isochrones==1.1.1':
	print('UPDATING '+starsdir+'/populations.py to remove reference to TESS band.')

	stars_popsfile = open(starsdir+'/populations.py', mode='r')
	new_stars_popsfile = open(starsdir+'/populations_new.py', mode='w')

	for nline, line in enumerate(stars_popsfile):
		#newline = line.replace('representation=', 'representation_type=')
		if line.startswith('BANDS'):
			newline = "BANDS = ['g','r','i','z','J','H','K','Kepler']\n" #### removed 'TESS' band, since it doesn't play with isochrones 1.1.1
		else:
			newline = line
		new_stars_popsfile.write(newline)

	stars_popsfile.close()
	new_stars_popsfile.close()
	os.system('mv '+starsdir+'/populations.py '+starsdir+'/populations_DEPRECATED.py')
	os.system('mv '+starsdir+'/populations_new.py '+starsdir+'/populations.py')

	print('UPDATED '+starsdir+'/populations.py')



#### UPDATE UTILS.PY
utilsfile = open(orbitsdir+'/utils.py', mode='r')
new_utilsfile = open(orbitsdir+'/utils_new.py', mode='w')

for nline, line in enumerate(utilsfile):
	newline = line.replace('representation=', 'representation_type=')
	new_utilsfile.write(newline)

utilsfile.close()
new_utilsfile.close()
os.system('mv '+orbitsdir+'/utils.py '+orbitsdir+'/utils_DEPRECATED.py')
os.system('mv '+orbitsdir+'/utils_new.py '+orbitsdir+'/utils.py')

print('UPDATED '+orbitsdir+'/utils.py')
print(' ')

print(' ')
print('UPDATED VESPA SCRIPTS! ')
print('can be found at '+orbitsdir+'/utils.py and populations.py')