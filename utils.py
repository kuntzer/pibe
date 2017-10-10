import cPickle as pickle
import gzip
import os
import numpy as np
import csv

def writepickle(obj, filepath, protocol = -1):
	"""
	I write your python object obj into a pickle file at filepath.
	If filepath ends with .gz, I'll use gzip to compress the pickle.
	Leave protocol = -1 : I'll use the latest binary protocol of pickle.
	"""
	if os.path.splitext(filepath)[1] == ".gz":
		pkl_file = gzip.open(filepath, 'wb')
	else:
		pkl_file = open(filepath, 'wb')
	
	pickle.dump(obj, pkl_file, protocol)
	pkl_file.close()
	
def readpickle(filepath):
	"""
	I read a pickle file and return whatever object it contains.
	If the filepath ends with .gz, I'll unzip the pickle file.
	"""
	if os.path.splitext(filepath)[1] == ".gz":
		pkl_file = gzip.open(filepath,'rb')
	else:
		pkl_file = open(filepath, 'rb')
	obj = pickle.load(pkl_file)
	pkl_file.close()
	return obj

def colnorm(cat, name, oname=None):
	"""
	Normalises the column `name` of astropy table `cat`. Optionnally outputs the normalised values to column `oname`. 
	"""
	min_ = cat[name].min()
	max_ = cat[name].max()
	
	col = (cat[name] - min_) / (max_ - min_)
	
	if oname is None:
		cat[name] = col
	else:
		cat[oname] = col
		
def find_nearest(array,value):
	''' Find nearest value is an array '''
	idx = (np.abs(array-value)).argmin()
	return idx

def load_bmg(fname, main_sequence):	
	data=[]
	with open(fname+'.dat') as observability_file:
		observability_data = csv.reader(observability_file, delimiter="\t")
		for row in observability_data:
		# if line is empty, skip otherwise filter out the blank
			dline=row[0].split()
			if len(dline)==17 and not dline[6].isdigit():
				dline.insert(6, '0')
			if dline[0][0]=='#': continue
			data.append(np.asarray(dline, dtype=np.float))
			
	data=np.asarray(data)
	if main_sequence: 
		data=data[data[:,2] == 5] #Takes only main sequence stars
	
	return data
