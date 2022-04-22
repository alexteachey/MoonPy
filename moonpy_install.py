import sys
import subprocess
import time


if sys.platform == 'darwin':
	#### use moonpy_env_macOS.yml
	your_OS = 'macOS'
	environment_name = 'moonpy_env_macOS'
	install_command = 'conda env create --file moonpy_env_macOS.yml && python pathmaker.py'

elif (sys.platform == 'linux') or (sys.platform == 'linux2'):
	your_OS = 'linux'
	environment_name = 'moonpy_env_linux'
	install_command = 'conda env create --file moonpy_env_linux.yml && python pathmaker.py'
	

print('established you are running '+your_OS+'.')

print('You are about to create a new conda environment and install MoonPy.')
print('ENVIRONMENT NAME: ', environment_name)
print(' ')

continue_install = input('Do you want to continue? y/n: ')
if continue_install == 'y':
	subprocess.Popen(install_command).wait()
else:
	raise Exception('you opted not to continue installation.')


#### print all conda environments
subprocess.Popen('conda env list')
print(' ')
print(' ')
print('You have created a new environment for MoonPy called '+environment_name)
print("To use MoonPy, remember to type 'conda activate "+environment_name+"'. ")
print(' ')

activate_env_now = input('Do you want to activate this environment now? y/n: ')

if activate_env_now == 'y':
	subprocess.Popen('conda activate '+environment_name)


