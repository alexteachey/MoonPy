import os


print('This script will delete an unwanted conda directory.')
print(' ')
os.system('conda env list')
print(' ')
print(' ')
env_to_remove = input('What is the name of the conda environment you want to remove? ')
os.system('conda env remove -n '+env_to_remove)
print('The conda environment '+env_to_remove+' has been removed.')