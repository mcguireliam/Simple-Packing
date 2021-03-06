'''
    the example of input is:
    python pdb2mrc.py ./pdbfile ./mrcfile

	This function will use eman2 package and call the following command:
	e2pdb2mrc.py ./pdbfile/3hhb.pdb ./mrcfile/3hhb.mrc res=2.8 apix=1.0
'''


import os
import sys

def pdb2mrc(pdbpath='../IOfile/pdbfile',mrcpath='../IOfile/mrcfile',res=10.0,apix=3.0):
	'''

	read all the pdb files in pdbpath and conver them into mrc files.
	the output files will be saved in mrcpath

	:param pdbpath: the path of input file.
	:param mrcpath: the path of output file
	:param res: the resolution of density map. Resolution defined in 'EM' way. reciprocal of the 1/2 width of a Gaussian in Fourier space
	:param apix: Angstroms per voxel in the output, default is 1
	:return: the mrc file of all the pdf file in pdbpath will save in mrcpath
	'''

	print('Using apix='+str(apix)+' and res='+str(res))
	# os.chdir(pdbpath)
	for filename in os.listdir(pdbpath):
		print(filename)
		if filename.endswith(".pdb"):
			os.system('e2pdb2mrc.py '+ pdbpath + '/' + filename + ' ' + mrcpath + '/' + filename[:-4]+'.mrc res='+str(res)+ ' apix='+str(apix))
	if any(".pdb" in s for s in os.listdir(pdbpath)) is False:
		print("ERROR: No pdb files found in path")
		return ValueError
	return None

if __name__ == '__main__':
	try:
		pdb2mrc(sys.argv[1],sys.argv[2])
	except:
		pdb2mrc()
