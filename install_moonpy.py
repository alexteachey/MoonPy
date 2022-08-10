import sys
import subprocess
import time
import os
from pathlib import Path 
import traceback
import time

homepath = str(Path.home())

possible_paths = sys.path
site_packages_paths = []
for possible_path in possible_paths:
	if ('site-packages' in possible_path) and ('conda' in possible_path) and ('site-packages/' not in possible_path):
		site_packages_path = possible_path 
		site_packages_paths.append(site_packages_path)
		print('site-packages path: ', site_packages_path)
if len(site_packages_paths) > 1:
	print("WARNING: multiple site-packages paths found. Pandora installation may be affected.")
	time.sleep(3)

def build_env_and_install(packagename, standard_environment_name):

	standard_environment_yml = 'env_setup_files/'+standard_environment_name+'.yml'

	print('You are about to create a new conda environment and install '+packagename)

	user_environment_name = input('What do you want to call the environment for '+packagename+'? [press ENTER to use '+standard_environment_name+']: ')
	if len(user_environment_name) == 0:
		user_environment_name = standard_environment_name 
	user_environment_yml = 'env_setup_files/'+user_environment_name+'.yml'

	if user_environment_name != standard_environment_name:
		print("Using USER environment name: ", user_environment_name)
		#### now we're going to alter that file
		env_file = open(standard_environment_yml, mode='r')
		new_env_file = open(user_environment_yml, mode='w')

		for nline,line in enumerate(env_file):
			#if nline == 0:
			if line.startswith('name') or line.startswith('prefix'):
				#newline = 'name: '+user_environment_name+'\n'
				newline = line.replace(standard_environment_name, user_environment_name)
			else:
				newline = line
			#newline = line.replace(standard_environment_name, user_environment_name)
			new_env_file.write(newline)

		env_file.close()
		new_env_file.close()

		environment_name = user_environment_name
		environment_yml = 'env_setup_files/'+user_environment_name+'.yml'

		print('environment_name: ', environment_name)
		print('environment_yml: ', environment_yml)

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
	if your_OS == 'linux':
		print(' ')
		print(' ')
		print("LINUX USERS: first 'conda activate "+environment_name+"', then type 'python pathmaker.py to complete installation.'")
		if packagename.lower() == 'vespa':
			print("IMPORTANT: before using VESPA, you need to update some scripts with deprecated keywords.")
			print("In the terminal, within the MoonPy directory and in the conda environment you've created for VESPA,")
			print("type 'python vespa_script_updater.py'. Do this now please.")
		print(' ')
		print(' ')

	return environment_name 


print(' ')
print(' ')
print(' ')



if sys.platform == 'darwin':
	#### use moonpy_env_macOS.yml
	your_OS = 'macOS'

	install_pandora_and_ultranest = input("Do you want to install Pandora and Ultranest? y/n: ")

	if install_pandora_and_ultranest == 'y': 
		standard_moonpy_environment_name = 'moonpy_with_pandora_macOS'
		backup_moonpy_environment_name = 'moonpy_env_macOS'
	else:
		standard_moonpy_environment_name = 'moonpy_env_macOS'
		backup_moonpy_environment_name = 'moonpy_env_macOS' ### same -- there is no backup

	standard_vespa_env_name = 'vespa_for_mac'
	#install_command = 'conda env create --file moonpy_env_macOS.yml'


elif (sys.platform == 'linux') or (sys.platform == 'linux2'):
	your_OS = 'linux'
	#standard_moonpy_environment_name = 'moonpy_env_linux'
	standard_moonpy_environment_name = 'moonpy_env_linux'
	backup_moonpy_environment_name = 'OLD_moonpy_env_linux' 
	standard_vespa_env_name = 'vespa_for_linux'
	#install_command = 'conda env create --file moonpy_env_linux.yml'

standard_moonpy_environment_yml = 'env_setup_files/'+standard_moonpy_environment_name+'.yml'
backup_moonpy_environment_yml = 'env_setup_files/'+backup_moonpy_environment_name+'.yml'
standard_vespa_env_yml = 'env_setup_files/'+standard_vespa_env_name+'.yml'

	
print(' ')
print(' ')
print('established you are running '+your_OS+'.')




### INSTALL MOONPY
try:
	moonpy_envname = build_env_and_install(packagename='MoonPy', standard_environment_name=standard_moonpy_environment_name)
except:
	moonpy_envname = build_env_and_install(packagename='MoonPy', standard_environment_name=backup_moonpy_environment_name)




#### VESPA INSTALLATION 
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
	vespa_envname = build_env_and_install(packagename='vespa', standard_environment_name=standard_vespa_env_name)








