import numpy as np 
from sklearn import decomposition
from astropy.io import fits
from astropy.table import Table, vstack
import os 

import utils

###################################################################################################
#

n = range(1)

in_dir = "output/datasets_packaged"

run_on = "train"

run_on_snr = "no"

do_pca = True

###################################################################################################
#

loaddir = os.path.join(in_dir, run_on, "snr_{}".format(run_on_snr))

training_source = None

datastamps = []
for ii in range(1):
	training_data_tmp = fits.getdata(os.path.join(loaddir, "{}_{:03}.fits".format(run_on, ii)))
	training_source_cat = Table.read(os.path.join(in_dir, run_on, "catalogs", "{}_{:03}_truth_cat.fits".format(run_on, ii)))
	
	counts = 0
	for xcol, ycol in training_source_cat["xpycat", "ypycat"]:
		stamp = training_data_tmp[xcol - 24: xcol + 24, ycol - 24: ycol + 24].flatten()
		datastamps.append(stamp)
		counts += 1
	print counts
		
	if training_source is None:
		training_source = training_source_cat
	else:
		training_source = vstack([training_source, training_source_cat])

datastamps = np.array(datastamps)
if do_pca:
	pca = decomposition.PCA(n_components=8, whiten=True)
	pca.fit(datastamps)
	utils.writepickle(pca, "pca.pkl")
else:
	pca = utils.readpickle("pca.pkl")

X = pca.transform(datastamps)
print X.shape

hh = X[:,0]
print hh.shape 
print np.unique(hh).shape
#hh = training_source["fwhm"]

import pylab as plt
plt.figure()
plt.hist(hh, 1000)

plt.figure()
plt.scatter(training_source["xfield"], training_source["yfield"], c=hh, edgecolor="None")
plt.colorbar()

plt.show()


