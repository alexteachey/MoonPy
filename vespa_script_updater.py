import sys 
import os


#### FIX VESPA -- OUT OF DATE KEYWORDS FROM ASTROPY
### '/opt/miniconda3/envs/vespa_with_moonpy_for_mac/lib/python3.6/site-packages/vespa/orbits'
pythonpaths = sys.path 

for pythonpath in pythonpaths:
	if 'site-packages' in pythonpath:
		desired_path = pythonpath

orbitsdir = desired_path+'/vespa/orbits'

#### UPDATE POPULATIONS.PY
popsfile = open(orbitsdir+'/populations.py', mode='r')
new_popsfile = open(orbitsdir+'/populations_new.py', mode='w')

for nline, line in enumerate(popsfile):
	newline = line.replace('representation=', 'representation_type=')
	new_popsfile.write(newline)

popsfile.close()
new_popsfile.close()
os.system('mv '+orbitsdir+'/populations.py '+orbitsdir+'/populations_DEPRECATED.py')
os.system('mv '+orbitsdir+'/populations_new.py '+orbitsdir+'/populations.py')



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

print(' ')
print('UPDATED VESPA SCRIPTS! ')
print('can be found at '+orbitsdir+'/utils.py and populations.py')