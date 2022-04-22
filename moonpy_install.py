import sys
import subprocess


if sys.platform == 'darwin':
	#### use moonpy_env_macOS.yml
	environment_name = 'moonpy_env_macOS'
	install_command = 'conda env create --file moonpy_env_macOS.yml && python pathmaker.py'

elif (sys.platform == 'linux') or (sys.platform == 'linux2'):
	environment_name = 'moonpy_env_linux'
	install_command = 'conda env create --file moonpy_env_linux.yml && python pathmaker.py'
	
subprocess.Popen(install_command).wait()


#### print all conda environments
subprocess.Popen('conda env list')
print(' ')
print(' ')
print('You have created a new environment for MoonPy called '+environment_name)

activate_env_now = input('Do you want to activate this environment now? y/n: ')

if activate_env_now == 'y':
	subprocess.Popen('conda activate '+environment_name)


