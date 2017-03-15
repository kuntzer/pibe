import cPickle as pickle
import gzip
import os

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
