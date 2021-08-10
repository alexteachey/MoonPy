import sys
import os

#moonpydir = os.path.realpath(__file__)
#moonpydir = moonpydir[:moonpydir.find('/pathmaker.py')]

def make_pathfile(moonpydir):
	#### create the moonpy.pth file
	os.system('touch '+moonpydir+'/moonpy.pth')
	moonpy_pthfile = open(moonpydir+'/moonpy.pth', mode='w')
	moonpy_pthfile.write(moonpydir)
	moonpy_pthfile.close()

	pythonpath = sys.executable
	prepath = pythonpath[:pythonpath.find('/bin')]
	libpath = prepath+'/lib'
	libfiles = os.listdir(libpath)
	for libfile in libfiles:
		if libfile.startswith('python'):
			pythonlib = libfile
			break
	pythonlibpath = libpath+'/'+pythonlib 
	site_packages_path = pythonlibpath+'/site-packages'
	source_path = moonpydir+'/moonpy.pth'
	destination_path = site_packages_path+'/moonpy.pth'

	#### copy moonpy.pth to the site_packages_path
	if os.path.exists(destination_path):
		print(destination_path+' exists.')
	else:
		os.system('cp '+source_path+' '+destination_path)
		print('copied '+source_path+' to '+destination_path)

make_pathfile(os.getcwd())