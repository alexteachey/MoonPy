import os

### identify your working environment
current_environment = os.environ['CONDA_DEFAULT_ENV']

### export this environment to a .yml file

output_name = input("What name do you want to give the output environment file? Press ENTER to keep it the same as the current environment.")

if output_name == '':
	output_name = current_environment

export_command = 'conda env export -n '+current_environment+' -f '+output_name+'.yml --no-builds'

os.system(export_command)

