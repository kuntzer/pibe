import galsim
import pylab as plt
import numpy as np
import os
import astropy.io.fits as fits
import datetime
import multiprocessing

import utils

###################################################################################################
# Defining variables
# Pixel scale in arcsec / pixel
pixel_scale = 0.1#/12.

image_size = 48

n_ini = 0
n_psf = 250# 20000#0

# Where to save the PSFs
out_dir = 'output/psf_nonoise'

# Where are the interpolation functions
interp_dir = "inputdata/fullpsffield"

# Directory with the spectra (where the README of Pickles is)
spectra_dir = "inputdata/spectra/pickles"

# 'STEP_04750x09700', # Step euclid in Angstrom
filterband_min = 4750. 
filterband_max = 9700.

# Which spectra for the star?
spectra_id = 27 # G5V

# Skip already drawn images in the output dir?
skipdone = True

# Show for debug purposes
show = False

# ncpu
ncpu = 5

if not os.path.exists(out_dir):
	os.makedirs(out_dir)

# GalSim paramters
psf_obsc = 0.29 # (0.35m / 1.2m) = 0.291667 # This gies a FWHM of 1" if not corrected
psf_nstruts = 3
psf_strut_thick = 0.015
psf_strut_angle = (360-105) * galsim.degrees

psf_defocus = .0	 # The aberrations are all taken to be quite modest here.
psf_astig1 = 0.065 * 0.8	  # (I don't actually know what are appropriate for HST...)
psf_astig2 = -0.08 * 0.8#6
psf_coma1 = 0.#2
psf_coma2 = 0.
psf_trefoil1 = -0.04 * 0.8#5#1
psf_trefoil2 = +0.01 * 0.8#1

###################################################################################################
# Initialisation
spectrum = fits.getdata(os.path.join(spectra_dir, "dat_uvi/pickles_%d.fits" % spectra_id), view=np.ndarray)
spectrum = spectrum.astype([('WAVELENGTH', '>f4'), ('FLUX', '>f4')]).view('>f4').reshape(len(spectrum), -1)

ids = np.where(spectrum[:,0] >= filterband_min)[0]
ids2 = np.where(spectrum[ids,0] <= filterband_max)[0]
ids = ids[ids2][::40]

if not os.path.exists("output"):
	os.mkdir("output")
	
ing1, ing2, infwhm = utils.readpickle(os.path.join(interp_dir, 'interp.pkl'))

a=np.linspace(0, 1, 100)
b=np.linspace(0, 1, 100)
c,d = np.meshgrid(a,b)
maxr = infwhm.ev(c,d).max()
minr = infwhm.ev(c,d).min()
maxg1 = ing1.ev(c,d).max()
ming1 = ing1.ev(c,d).min()
maxg2 = ing2.ev(c,d).max()
ming2 = ing2.ev(c,d).min()

print 'Euclid limits'
print 'g1', ming1, maxg1
print 'g2', ming2, maxg2
print 'fwhm', minr, maxr
print '****************'

def corr_size(r, minr, maxr, effmaxr=.3):
	return (r - minr) / (maxr - minr) * effmaxr + 1.

###################################################################################################
# defining the worker function
def worker(params):
	
	# Making sure that the seed isn't the same for all processes
	np.random.seed()
	
	istar = params
	
	imfn = os.path.join(out_dir, "star_%05d.fits" % istar)
	
	if os.path.exists(imfn) and skipdone:
		print 'star {} already exists, skipping'.format(istar)
		return

	xstar, ystar = np.random.uniform(size=(2))
	
	g1 = ing1(xstar, ystar)[0][0]
	g2 = ing2(xstar, ystar)[0][0]
	fwhm = infwhm(xstar, ystar)[0][0]
	
	print 'Star {}: (x={:.3f},y={:.3f})'.format(istar, xstar, ystar)
	psf_imgi = None
	then = datetime.datetime.now()
	for ii, ilam in enumerate(ids):
	
		lam = spectrum[ilam,0] # unit: angstrom
		flux = spectrum[ilam,1] # unit: flam
		
		if show: print '\tWavelength ({})'.format(istar), ii+1, '/', len(ids), " : ", lam*0.1, "nm"
		
		psf = galsim.OpticalPSF(
	#		lam_over_diam=psf_lam_over_D, 
			lam=lam*0.1, diam=1.2, # Galsim works in nm, diam in m
			obscuration=psf_obsc,# * corr_size(fwhm, minr, maxr),
			nstruts=psf_nstruts, strut_thick=psf_strut_thick, strut_angle=psf_strut_angle,
			defocus=psf_defocus, astig1=psf_astig1* corr_size(fwhm, minr, maxr), astig2=psf_astig2* corr_size(fwhm, minr, maxr),
			coma1=psf_coma1* corr_size(fwhm, minr, maxr), coma2=psf_coma2* corr_size(fwhm, minr, maxr), trefoil1=psf_trefoil1* corr_size(g1, ming1, maxg1),
			 trefoil2=psf_trefoil2* corr_size(g1, ming1, maxg1))

		#print '%d\t%+1.2f\t%+1.2f' % (i, g1, g2)
		
		psf = psf.shear(g1=g1/2., g2=g2*1.7)
		
		image = galsim.ImageF(image_size, image_size)
		psf.drawImage(image=image, scale=pixel_scale)

		#if image.array.min() < 0:
		#	print "WARNING: <0 values!, min is", image.array.min()
		#	image.array[image.array < 0] = 0 # positivity constraint
			
		slicei = image.array * flux 
		if psf_imgi is None:
			psf_imgi = slicei
		else:
			psf_imgi += slicei
			
	psf_imgi /= np.sum(spectrum[ids,1])
	
	now = datetime.datetime.now()
	
	h = fits.Header()
	h['x'] = xstar
	h['y'] = ystar
	
	#TODO: jitter here!
	#if jiiter:
	#	ud = galsim.UniformDeviate() # This gives a random float in [0, 1)
	#	# We apply some jitter to the position of this psf
	#	xjitter = ud() - 0.5 # This is the minimum amount -- should we do more, as real stars are not that well centered in their stamps ?
	#	yjitter = ud() - 0.5
	#	psf_imgi = psf_imgi.shift(xjitter,yjitter)
	
	fits.writeto(imfn, psf_imgi, clobber=True, header=h)

	print 'Star {}: done, took {:s}'.format(istar, now - then)
	if show:
		plt.figure()
		plt.imshow(np.log10(psf_imgi), interpolation='None')
		plt.title(istar)
		plt.show()
	
	
###################################################################################################
# Running the worker
params = range(n_ini, n_psf)

if ncpu == 1: # The single-processing version, much easier to debug !
	makeimg = map(worker, params)

else: # The simple multiprocessing map is:
	pool = multiprocessing.Pool(processes=ncpu)
	makeimg = pool.map(worker, params)
	pool.close()
	pool.join()
		

