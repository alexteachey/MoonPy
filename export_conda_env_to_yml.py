import os

### identify your working environment
current_environment = os.environ['CONDA_DEFAULT_ENV']

### export this environment to a .yml file

output_name = input("What name do you want to give the output environment file? (Press ENTER to keep it the same as the current environment): ")

if output_name == '':
	output_name = current_environment

export_command = 'conda env export -n '+current_environment+' -f '+output_name+'.yml --no-builds'

os.system(export_command)



#### now make sure you replace the current_environment name with the output name!
exported_filename = output_name+'.yml'
replacement_filename = 'replacement_'+output_name+'.yml'

orig_file = open(exported_filename, mode='r')
new_file = open(replacement_filename, mode='w')

for nline,line in enumerate(orig_file):
	line = line.replace(current_environment, output_name)
	new_file.write(line)

orig_file.close()
new_file.close()

### replace the original file with the new file
os.system('mv '+replacement_filename+' '+exported_filename)
