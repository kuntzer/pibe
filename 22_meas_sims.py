from astropy.table import Table
from astropy.io import fits
import galsim
import glob
import numpy as np
import os
import pylab as plt

###################################################################################################
# Defining variables

# output dir
sim_dir = "output/psf_nonoise"

###################################################################################################
# Initialisation

psfs = glob.glob(os.path.join(sim_dir, "star_*.fits"))

print "Found {} simulations".format(len(psfs))

og1 = []
og2 = []
or2 = []
ox = []
oy = []

###################################################################################################
# measurements of the PSFs

for istar, fnpsf in enumerate(psfs):
	g5v_img = None
	
	fit = fits.open(fnpsf)[0]

	data = fit.data
	header = fit.header

	# Moving all of this to GalSim format for measurement
	data = galsim.ImageF(data)
	res = galsim.hsm.FindAdaptiveMom(data)
	
	ox.append(header['X'])
	oy.append(header['Y'])
	og1.append(res.observed_shape.g1)
	og2.append(res.observed_shape.g2)
	or2.append(res.moments_sigma * np.sqrt(2. *np.log(2)) * 2. / 12.)
	
	print "meas:g1={:+.4f}\tg2={:+.4f}\tFWHM={:+.4f}".format( \
		res.observed_shape.g1, res.observed_shape.g2, res.moments_sigma * np.sqrt(2. *np.log(2)) * 2. / 12.)
	
	print "true:g1={:+.4f}\tg2={:+.4f}\tFWHM={:+.4f}".format(header['G1'], header['G2'], header["FWHM"])
	print
	
plt.figure()
cb = plt.scatter(ox, oy, c=og1, edgecolors="None")
plt.colorbar(cb)
plt.title("g1")

plt.figure()
cb = plt.scatter(ox, oy, c=og2, edgecolors="None")
plt.colorbar(cb)
plt.title("g2")

plt.figure()
cb = plt.scatter(ox, oy, c=or2, edgecolors="None")
plt.colorbar(cb)
plt.title("fwhm")

plt.show()
