import numpy as np 
from sklearn import decomposition
from astropy.io import fits
from astropy.table import Table, vstack
import os 
import itertools

import utils

###################################################################################################
#

n = range(3)

in_dir = "output/datasets_packaged"
out_dir = "output/encode_pca"

run_on = "train"

run_on_snr = "no"

do_pca = True

###################################################################################################
#

loaddir = os.path.join(in_dir, run_on, "snr_{}".format(run_on_snr))

training_source = None

datastamps = []
for ii in n:
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

datasets = ["train", "test"]
snrs = ["no", "20", "100"]
zzip = list(itertools.product(*[datasets, snrs]))

if not os.path.exists(out_dir):
	os.mkdir(out_dir)

for ro, rs in zzip:
	
	loaddir = os.path.join(in_dir, ro, "snr_{}".format(rs))
	training_source = None
	
	print "Working in:", loaddir
	
	if ro == "train":
		nps = 5
	elif ro == "test":
		nps = 10
		
	for ii in range(nps):
		datastamps = []
		
		data = fits.getdata(os.path.join(loaddir, "{}_{:03}.fits".format(ro, ii)))
		source_cat = Table.read(os.path.join(in_dir, ro, "catalogs", "{}_{:03}_truth_cat.fits".format(ro, ii)))
		
		counts = 0
		for xcol, ycol in source_cat["xpycat", "ypycat"]:
			stamp = data[xcol - 24: xcol + 24, ycol - 24: ycol + 24].flatten()
			datastamps.append(stamp)
			counts += 1
		print counts
			
		datastamps = np.array(datastamps)

		codes = pca.transform(datastamps)
		xys = np.vstack([np.array(source_cat["xfield"]), np.array(source_cat["yfield"])]).T
		codes = np.hstack([xys, codes])
		
		utils.writepickle(codes, os.path.join(out_dir, "snr_{}_{}_{:03}.pkl".format(rs, ro, ii)))
		

