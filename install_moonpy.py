import sys
import subprocess
import time
import os


if sys.platform == 'darwin':
	#### use moonpy_env_macOS.yml
	your_OS = 'macOS'
	standard_environment_name = 'moonpy_env_macOS'
	#install_command = 'conda env create --file moonpy_env_macOS.yml'

elif (sys.platform == 'linux') or (sys.platform == 'linux2'):
	your_OS = 'linux'
	standard_environment_name = 'moonpy_env_linux'
	#install_command = 'conda env create --file moonpy_env_linux.yml'
standard_environment_yml = standard_environment_name+'.yml'

	
print(' ')
print(' ')
print('established you are running '+your_OS+'.')

print('You are about to create a new conda environment and install MoonPy.')

user_environment_name = input('What do you want to call the environment? [press ENTER to use '+standard_environment_name+']: ')
if len(user_environment_name) == 0:
	user_environment_name = standard_environment_name 
user_environment_yml = user_environment_name+'.yml'

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
print('creating the MoonPy path...')
subprocess.run('source activate '+environment_name+' && python pathmaker.py && conda deactivate', shell=True)


#### print all conda environments
try:
	subprocess.Popen('conda env list', shell=True).wait()
except:
	try:
		subprocess.Popen('conda env list', shell=False).wait()
	except:
		os.system('conda env list')
print(' ')
print(' ')
print('You have created a new environment for MoonPy called '+environment_name)
print("To use MoonPy, remember to type 'conda activate "+environment_name+"'. ")
print('You can than import moonpy like any other python package.')






