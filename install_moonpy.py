import sys
import subprocess
import time
import os
from pathlib import Path 
import traceback

homepath = str(Path.home())

possible_paths = sys.path
site_packages_paths = []
for possible_path in possible_paths:
	if 'site-packages' in possible_path:
		site_packages_path = possible_path 
		site_packages_paths.append(site_packages_path)
		print('site-packages path: ', site_packages_path)
if len(site_packages_paths) > 1:
	print("WARNING: multiple site-packages paths found. Pandora installation may be affected.")


def build_env_and_install(packagename, standard_environment_name):
	standard_environment_yml = 'env_setup_files/'+standard_environment_name+'.yml'

	print('You are about to create a new conda environment and install '+packagename)

	user_environment_name = input('What do you want to call the environment for '+packagename+'? [press ENTER to use '+standard_environment_name+']: ')
	if len(user_environment_name) == 0:
		user_environment_name = standard_environment_name 
	user_environment_yml = 'env_setup_files/'+user_environment_name+'.yml'

	if user_environment_name != standard_environment_name:
		#### now we're going to alter that file
		env_file = open(standard_environment_yml, mode='r')
		new_env_file = open(user_environment_yml, mode='w')

		for nline,line in enumerate(env_file):
			newline = line.replace(standard_environment_name, user_environment_name)
			new_env_file.write(newline)

		env_file.close()
		new_env_file.close()

		environment_name = user_environment_name
		environment_yml = user_environment_yml

	else:
		environment_name = standard_environment_name
		environment_yml = standard_environment_yml 

	install_command = 'conda env create --file '+environment_yml


	print('ENVIRONMENT NAME: ', environment_name)
	print(' ')

	continue_install = input('Do you want to continue? y/n: ')
	if continue_install == 'y':
		subprocess.Popen(install_command, shell=True).wait()
	else:
		raise Exception('you opted not to continue installation.')


	#### activate environment, build the moonpy path, and deactivate 
	print('creating the '+packagename+' path...')
	if packagename.lower() == 'moonpy':
		#### only need to do this for moonpy.
		result = subprocess.run('source activate '+environment_name+' && python pathmaker.py && conda deactivate', shell=True, capture_output=True, text=True)


	elif packagename.lower() == 'vespa':
		#### run the vespa updater 
		result = subprocess.run('source activate '+environment_name+' && python pathmaker.py && python vespa_script_updater.py && conda deactivate', shell=True, capture_output=True, text=True)

	print('stdout: ', result.stdout)
	print('stderr: ', result.stderr)		

	try:
		subprocess.Popen('conda env list', shell=True).wait()
	except:
		try:
			subprocess.Popen('conda env list', shell=False).wait()
		except:
			os.system('conda env list')
	print(' ')
	print(' ')
	print('You have created a new environment for '+packagename+' called '+environment_name)
	print("To use "+packagename+", remember to type 'conda activate "+environment_name+"'. ")
	print('You can now import it like any other python package.')
print(' ')
print(' ')
print(' ')



if sys.platform == 'darwin':
	#### use moonpy_env_macOS.yml
	your_OS = 'macOS'
	standard_moonpy_environment_name = 'moonpy_env_macOS'
	standard_vespa_env_name = 'vespa_for_mac'
	#install_command = 'conda env create --file moonpy_env_macOS.yml'

elif (sys.platform == 'linux') or (sys.platform == 'linux2'):
	your_OS = 'linux'
	standard_moonpy_environment_name = 'moonpy_env_linux'
	standard_vespa_env_name = 'vespa_for_linux'
	#install_command = 'conda env create --file moonpy_env_linux.yml'
standard_moonpy_environment_yml = 'env_setup_files/'+standard_moonpy_environment_name+'.yml'
standard_vespa_env_yml = 'env_setup_files/'+standard_vespa_env_name+'.yml'

	
print(' ')
print(' ')
print('established you are running '+your_OS+'.')




### INSTALL MOONPY
build_env_and_install(packagename='MoonPy', standard_environment_name=standard_moonpy_environment_name)


#### set up Pandora
print('Attempting to install up PANDORA...')
try:
	if your_OS == 'macOS':
		#### should just be able to pip install directly
		subprocess.run('pip install pandoramoon', shell=True, capture_output=True, text=True)

	elif your_OS == 'linux':
		#### need to delete a file first
		for site_packages_path in site_packages_paths:
			#### remove this conflicting file, and then pip install
			subprocess.run('rm -f '+site_packages_path+'/llvmlite*egg-info && pip install pandoramoon', shell=True, capture_output=True, text=True)
except:
	traceback.print_exc()
	print(' ')
	print(' ')
	print('PANDORA installation failed.')
	print(' ')
	print(' ')











setup_vespa = input("Do you want to install Tim Morton's VESPA code (recommended)? y/n: ")

if (setup_vespa == 'y') or (setup_vespa == ''):

	if os.path.exists(homepath+'/.isochrones'):
		print(' ')
	
		print('It is recommended that any existing version of ~/.isochrones be clobbered to ensure proper performance.')
		print('We propose instead to rename it in case of an installation error. ')
		print(' ')
		rename_or_clobber_or_nothing = input("Do you want to 'r'ename, 'c'lobber, or do 'n'othing with ~/.isochrones? ")
		if rename_or_clobber_or_nothing.lower() == 'r':
			os.system('mv '+homepath+'/.isochrones '+homepath+'/.isochrones_BACKUP')
			print('renamed ~/.isochrones to ~/.isochrones_BACKUP')
			print('isochrones will be downloaded automatically upon first run of VESPA.')

		elif rename_or_clobber_or_nothing.lower() == 'c':
			os.system('rm -rf '+homepath+'/.isochrones')
			print('removed ~/.isochrones')
			print('isochrones will be downloaded automatically upon first run of VESPA.')		
				
		elif rename_or_clobber_or_nothing.lower() == 'n':
			print('You opted to leave ~/.isochrones alone. If you run into problems, consider removing.')
		
		
	time.sleep(3)

	#### NOW INSTALL VESPA!
	build_env_and_install(packagename='vespa', standard_environment_name=standard_vespa_env_name)








