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
sim_dir = "output/psf_nonoise_smallpx_SED44"

###################################################################################################
# Initialisation

psfs = glob.glob(os.path.join(sim_dir, "star_*.fits"))

print "Found {} simulations".format(len(psfs))

g1s = []
g2s = []
fwhms = []
xstars = []
ystars = []
catseqs = []
colours = []

###################################################################################################
# measurements of the PSFs

for istar, fnpsf in enumerate(psfs):
	g5v_img = None
	
	catseqs.append(istar)
	
	fit = fits.open(fnpsf)[0]

	data = fit.data
	header = fit.header

	# Moving all of this to GalSim format for measurement
	data = galsim.ImageF(data)
	res = galsim.hsm.FindAdaptiveMom(data)
	
	xstars.append(header['X'])
	ystars.append(header['Y'])
	colours.append(header['spectrid'])
	g1s.append(res.observed_shape.g1)
	g2s.append(res.observed_shape.g2)
	fwhms.append(res.moments_sigma * np.sqrt(2. *np.log(2)) * 2. / 12.)
	
	print "meas:g1={:+.4f}\tg2={:+.4f}\tFWHM={:+.4f}".format( \
		res.observed_shape.g1, res.observed_shape.g2, res.moments_sigma * np.sqrt(2. *np.log(2)) * 2. / 12.)
	
plt.figure()
cb = plt.scatter(xstars, ystars, c=g1s, edgecolors="None")
plt.colorbar(cb)
plt.title("g1")

plt.figure()
cb = plt.scatter(xstars, ystars, c=g2s, edgecolors="None")
plt.colorbar(cb)
plt.title("g2")

plt.figure()
cb = plt.scatter(xstars, ystars, c=fwhms, edgecolors="None")
plt.colorbar(cb)
plt.title("fwhm")

	
cat = Table([catseqs, xstars, ystars, g1s, g2s, fwhms], names=('catseq', 'x', 'y', 'g1', 'g2', 'fwhm'))
cat['fwhm'].unit = 'arcsec'

catfn = os.path.join(sim_dir, "catsimmeas.fits")
cat.write(catfn, overwrite=True)

print cat

print "Catalog saved to {}".format(catfn)

plt.show()
