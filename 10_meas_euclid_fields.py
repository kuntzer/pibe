from astropy.table import Table
from astropy.io import fits
import galsim
import glob
import numpy as np
import os

###################################################################################################
# Defining variables

# 'STEP_05500x09000', # Step euclid in Angstrom
filterband_min = 5500. 
filterband_max = 6000.#9000.

# Which spectra should be used?
spectra_id = 27 # G5V

# Directory with the Euclid psf
psf_dir = "inputdata/inputpsf"

# Directory with the spectra (where the README of Pickles is)
spectra_dir = "inputdata/spectra/pickles"

# output dir
out_dir = "inputdata/psffield"

# Show star image before measurement? (only for debug purposes)
show = False

###################################################################################################
# Initialisation

euclid_psfs = glob.glob(os.path.join(psf_dir, "PSF_MC_T0072.ZMX__Npup_4096x4096_lbda_0.800um_central4ccds.cat_*.fits"))

print "Found {} Euclid PSFs files".format(len(euclid_psfs))

spectrum = fits.getdata(os.path.join(spectra_dir, 'dat_uvi/', "pickles_%d.fits" % spectra_id), view=np.ndarray)
spectrum = spectrum.astype([('WAVELENGTH', '>f4'), ('FLUX', '>f4')]).view('>f4').reshape(len(spectrum), -1)

total_flux_g5v = 0.
nslice = 0

xstars = []
ystars = []
wavls = []
catseqs = []
g1s = []
g2s = []
fwhms = []

if not os.path.exists(out_dir):
	os.mkdir(out_dir)
	
if show:
	import pylab as plt

###################################################################################################
# Generating the right spectra + measurements of the PSFs

for istar, fnpsf in enumerate(euclid_psfs):
	print 'Measuring star nb {}/{}'.format(istar+1, len(euclid_psfs))
	g5v_img = None
	
	cubepsf = fits.open(fnpsf)
	
	for ipsf, wpsf in enumerate(cubepsf):
		
		if ipsf == 0:
			continue
		elif ipsf == 1:
			xstars.append(wpsf.header["XFIELD"])
			ystars.append(wpsf.header["YFIELD"])
			wavl = wpsf.header["WLGTH0"] * 10000
			wavls.append(wavl)
			catseqs.append(wpsf.header["CATSEQ"])
			
		wavl = wpsf.header["WLGTH0"] * 10000
			
		if wavl < filterband_min or wavl > filterband_max :
			continue

		fspec = spectrum[:,0] == wavl
		# Check that we found one matching wavelenght
		assert np.size(spectrum[fspec, 1]) == 1
		fspec = spectrum[fspec, 1]
		
		if istar == 0:
			total_flux_g5v += fspec
			nslice += 1
			
		image_slice = wpsf.data * fspec
		if g5v_img is None:
			g5v_img = image_slice
		else:
			g5v_img += image_slice
	
	# Normalising	
	g5v_img /= float(nslice)
		
	if show:
		plt.figure()
		plt.imshow(np.log10(g5v_img), interpolation="None")
		plt.show()

	# Moving all of this to GalSim format for measurement
	g5v_img = galsim.ImageF(g5v_img)
	res = galsim.hsm.FindAdaptiveMom(g5v_img)
	
	g1s.append(res.observed_shape.g1)
	g2s.append(res.observed_shape.g2) 
	fwhms.append(res.moments_sigma * np.sqrt(2. *np.log(2)) * 2. / 12.) # We divide by 12 to get the value in Euclid px
	
cat = Table([catseqs, xstars, ystars, g1s, g2s, fwhms], names=('catseq', 'xfield', 'yfield', 'g1', 'g2', 'fwhm'))
cat['fwhm'].unit = 'arcsec'

catfn = os.path.join(out_dir, "catmeas.fits")
cat.write(catfn, overwrite=True)

print cat

print "Catalog saved to {}".format(catfn)
