import numpy as np 
from sklearn import decomposition
from astropy.io import fits
from astropy.table import Table, vstack
import os 
import itertools
import scipy.interpolate as interp
import pylab as plt
try:
	from pykrige.ok import OrdinaryKriging
	import pykrige.kriging_tools as kt
	krige = True
except:
	print "Kriging will not work: pykrige not available"
	krige = False
	
import utils

###################################################################################################
#

n = range(1)

in_dir = "output/datasets_packaged"
out_dir = "output/encode_pca"

run_on = "train"

run_on_snr = "20"

n_components = 8

do_pca = False
do_interp = True
save_interp = False
show = False

interpolantnames = ["SmoothBivariateSpline", "linRBF"]
if krige:
	interpolantnames.append("OrdinaryKriging")

###################################################################################################
# Learning to compress

if not os.path.exists(out_dir):
	os.mkdir(out_dir)

loaddir = os.path.join(in_dir, run_on, "snr_{}".format(run_on_snr))

training_source = None

datastamps = []
for ii in n:
	training_data_tmp = fits.getdata(os.path.join(loaddir, "{}_{:03}.fits".format(run_on, ii)))
	training_source_cat = Table.read(os.path.join(in_dir, run_on, "catalogs", "{}_{:03}_truth_cat.fits".format(run_on, ii)))
	
	counts = 0
	for xcol, ycol in training_source_cat["xpycat", "ypycat"]:
		stamp = training_data_tmp[int(xcol - 24): int(xcol + 24), int(ycol - 24): int(ycol + 24)].flatten()
		datastamps.append(stamp)
		counts += 1
	print counts
		
	if training_source is None:
		training_source = training_source_cat
	else:
		training_source = vstack([training_source, training_source_cat])

datastamps = np.array(datastamps)
if do_pca:
	print 'Fitting PCAs...'
	pca = decomposition.PCA(n_components=n_components, whiten=True)
	pca.fit(datastamps)
	utils.writepickle(pca, os.path.join(out_dir, "pca_{}.pkl".format(run_on_snr)))
	print 'Done.'
else:
	pca = utils.readpickle(os.path.join(out_dir, "pca_{}.pkl".format(run_on_snr)))

X = pca.transform(datastamps)

###################################################################################################
# Learning to fly... (No just kidding, learning to interpolate)

if do_interp:
	interpolants = {}
	for intername in interpolantnames:
		print "Learning interpolant **{}**".format(intername)
		
		x = training_source["xfield"]
		y = training_source["yfield"]
		if intername == "SmoothBivariateSpline":
			interpolants["SmoothBivariateSpline"] = []
			for ipca in range(n_components):

				print "Learning to interpolate component {}...".format(ipca)
				intw = interp.SmoothBivariateSpline(x, y, X[:,ipca])
				interpolants["SmoothBivariateSpline"].append(intw)
		elif intername == "OrdinaryKriging":
			interpolants["OrdinaryKriging"] = []
			for ipca in range(n_components):
				x = training_source["xfield"]
				y = training_source["yfield"]
				print "Learning to interpolate component {}...".format(ipca)
				intw = OrdinaryKriging(x, y, X[:,ipca], variogram_model='linear',
                     verbose=False, enable_plotting=False)
				interpolants["OrdinaryKriging"].append(intw)
		elif intername == "linRBF":
			interpolants["linRBF"] = []
			for ipca in range(n_components):
				x = training_source["xfield"]
				y = training_source["yfield"]
				print "Learning to interpolate component {}...".format(ipca)
				intw = interp.Rbf(x, y, X[:,ipca], function='linear')
				interpolants["linRBF"].append(intw)
		else:
			raise RuntimeError("Unknown interpolation scheme: {}".format(intername))
	if save_interp:
		utils.writepickle(interpolants, os.path.join(out_dir, "interppca.pkl"))
else:
	interpolants = utils.readpickle(os.path.join(out_dir, "interppca.pkl"))

if show:
	for testpcacomp in range(8):
		print testpcacomp
		hh = X[:,testpcacomp]
		#hh = training_source["fwhm"]
		
		plt.figure()
		ax = plt.subplot(1,2,1)
		plt.scatter(training_source["xfield"], training_source["yfield"], c=hh, edgecolor="None")
		
		ax = plt.subplot(1,2,2)
		XX, YY = np.meshgrid(np.linspace(0, 1), np.linspace(0, 1))
		plt.contourf(XX, YY, interpolants[testpcacomp].ev(XX, YY), 100)
		plt.title("PCA component {}".format(testpcacomp))
	plt.show()

###################################################################################################
# Predictions

datasets = ["train", "test"]
snrs = [run_on_snr]
zzip = list(itertools.product(*[datasets, snrs]))

for intername in interpolantnames:
	if not os.path.exists(os.path.join(out_dir, intername)):
		os.mkdir(os.path.join(out_dir, intername))
	
	print "Predicting with interpolant **{}**".format(intername)
	
	for ro, rs in zzip:
		
		loaddir = os.path.join(in_dir, ro, "snr_{}".format(rs))
		training_source = None
		
		print "Working in:", loaddir
		
		if ro == "train":
			nps = 1#5
		elif ro == "test":
			nps = 1#0
			
		for ii in range(nps):
			datastamps = []
			
			data = fits.getdata(os.path.join(loaddir, "{}_{:03}.fits".format(ro, ii)))
			source_cat = Table.read(os.path.join(in_dir, ro, "catalogs", "{}_{:03}_truth_cat.fits".format(ro, ii)))
			x = source_cat["xfield"]
			y = source_cat["yfield"]
			
			for xcol, ycol in source_cat["xpycat", "ypycat"]:
				stamp = data[int(xcol - 24): int(xcol + 24), int(ycol - 24): int(ycol + 24)].flatten()
				datastamps.append(stamp)
			datastamps = np.asarray(datastamps)
			
			composants = pca.transform(datastamps)
			xys = np.vstack([np.array(source_cat["xfield"]), np.array(source_cat["yfield"])]).T
			composants = np.hstack([xys, composants])
			utils.writepickle(composants, os.path.join(out_dir, intername, "snr_{}_{}_calib_{:03}.pkl".format(rs, ro, ii)))
					
			
			composants = []
			for inti in interpolants[intername]:
				if intername == "SmoothBivariateSpline" :
					composants.append(inti.ev(x, y))
				elif intername == "linRBF" :
					composants.append(inti(x, y))
				elif intername == "OrdinaryKriging":
					compo = inti.execute("points", x, y)
					compo = compo[0]
					composants.append(compo)
			composants = np.array(composants).T
	
			composants = np.hstack([xys, composants])
			
			
			utils.writepickle(composants, os.path.join(out_dir, intername, "snr_{}_{}_{:03}.pkl".format(rs, ro, ii)))
		

