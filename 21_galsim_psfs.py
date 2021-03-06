import galsim
import pylab as plt
import numpy as np
import os
import astropy.io.fits as fits
import datetime
import multiprocessing
import re 

import utils

###################################################################################################
# Defining variables
# Pixel scale in arcsec / pixel
pixel_scale = 0.1#/12.

image_size = 48

n_ini = 0
n_psf = 2000#00#15#0000

# Where to save the PSFs
out_dir = 'output/psf_nonoise_euclidlike'

# Where are the interpolation functions
interp_dir = "inputdata/fullpsffield"

# Directory with the spectra (where the README of Pickles is)
spectra_dir = "inputdata/spectra/pickles"

# 'STEP_05500x09000', # Step euclid in Angstrom
filterband_min = 5500. 
filterband_max = 9000.

# Which spectra for the star?
# G5V == 27
"""
# FILENAME      SPTYPE    EFF TEMP. (K)
pickles_uk_1    O5V     39810.7
pickles_uk_2    O9V     35481.4
pickles_uk_3    B0V     28183.8
pickles_uk_4    B1V     22387.2
pickles_uk_5    B3V     19054.6
pickles_uk_6    B6V   14125.4 # B5-7V, but changed to be read automatically
pickles_uk_7    B8V     11749.0
pickles_uk_9    A0V     9549.93
pickles_uk_10   A2V     8912.51
pickles_uk_11   A3V     8790.23
pickles_uk_12   A5V     8491.80
pickles_uk_14   F0V     7211.08
pickles_uk_15   F2V     6776.42
pickles_uk_16   F5V     6531.31
pickles_uk_20   F8V     6039.48
pickles_uk_23   G0V     5807.64
pickles_uk_26   G2V     5636.38
pickles_uk_27   G5V     5584.70
pickles_uk_30   G8V     5333.35
pickles_uk_31   K0V     5188.00
pickles_uk_33   K2V     4886.52
pickles_uk_36   K5V     4187.94
pickles_uk_37   K7V     3999.45
pickles_uk_38   M0V     3801.89
pickles_uk_40   M2V     3548.13
pickles_uk_43   M4V     3111.72
pickles_uk_44   M5V     2951.21
"""
spectra_ids = [1,2,3,4,5,6,7,9,10,11,12,14,15,16,20,23,26,27,30,31,33,36,37,38,40,43,44] 

# Distribution of spectrum (either "flat" or give path BMG path [with all text commented or removed])
# Should be a list
distrib_spectrum = ["inputdata/BGM/star_field_BGM_i_180_15_{}".format(fid) for fid in range(1,6)]

# Skip already drawn images in the output dir?
skipdone = True

# Show for debug purposes
show = False

# Include jitter?
jitter = True

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

if type(distrib_spectrum) == list:
	try:
		catstar = utils.readpickle(out_dir+"_catstar.pkl")
		print 'loading from file'
	except:
		catstar = None
		for dsfn in  distrib_spectrum:
			print 'loading BMG file: {}'.format(dsfn)
		
			starss = utils.load_bmg(dsfn, main_sequence=True)
			
			if catstar is None:
				catstar = starss
			else:
				catstar = np.vstack([catstar, starss])
		
		ids = np.arange(np.shape(catstar)[0])
		np.random.shuffle(ids)
		catstar = catstar[ids]
		utils.writepickle(catstar, out_dir+"_catstar.pkl")
			
# Use spectra/pickles/UVILIB 
# Builds the reference array
# Stellar types are 1=O, ..., M = 7
stellar_types = np.array(['', 'O', 'B', 'A', 'F', 'G', 'K', 'M', 'AGB', 'WD'])

st_fnames = np.genfromtxt("inputdata/spectra/pickles/filenames.dat", dtype=["S15", "S15", "f8"])
spectra_fnames = np.empty([len(st_fnames), 3])
for ii, (f, st, _) in enumerate(st_fnames) :
	sclass = np.where(st[0] == stellar_types)[0][0]
	ssubclass = st[1:-1]
	m = re.search('(?<=uk_)\w+', f)
	spectra_fnames[ii] = [sclass, ssubclass, m.group(0)]

###################################################################################################
# defining the worker function
def worker(params):
	
	# Making sure that the seed isn't the same for all processes
	np.random.seed()
	
	istar = params
	if distrib_spectrum == "flat":
		spectra_id = np.random.choice(spectra_ids)
	elif type(distrib_spectrum) == list:
		spec_name = catstar[istar,3]
		
		sub, let = np.modf(spec_name)
		spec_dispo = spectra_fnames[spectra_fnames[:,0]==let,:]
		idsub = utils.find_nearest(spec_dispo[:,1], sub * 10.)
		spectra_id = spec_dispo[idsub,2]

	spectrum = fits.getdata(os.path.join(spectra_dir, "dat_uvi/pickles_%d.fits" % spectra_id), view=np.ndarray)
	spectrum = spectrum.astype([('WAVELENGTH', '>f4'), ('FLUX', '>f4')]).view('>f4').reshape(len(spectrum), -1)
	
	ids = np.where(spectrum[:,0] >= filterband_min)[0]
	ids2 = np.where(spectrum[ids,0] <= filterband_max)[0]
	ids = ids[ids2][::40]
	
	imfn = os.path.join(out_dir, "star_%05d.fits" % istar)
	
	if os.path.exists(imfn) and skipdone:
		print 'star {} already exists, skipping'.format(istar)
		return

	xstar, ystar = np.random.uniform(size=(2))
	
	g1 = ing1(xstar, ystar)[0][0]
	g2 = ing2(xstar, ystar)[0][0]
	fwhm = infwhm(xstar, ystar)[0][0]
	
	print 'Star {}: (x={:.3f},y={:.3f}), spectrum id={}'.format(istar, xstar, ystar, spectra_id)
	psf_imgi = None
	then = datetime.datetime.now()
	
	if jitter:
		ud = galsim.UniformDeviate() # This gives a random float in [0, 1)
		# We apply some jitter to the position of this psf
		xjitter = ud() - 0.5 # This is the minimum amount -- should we do more, as real stars are not that well centered in their stamps ?
		yjitter = ud() - 0.5
		xjitter *= pixel_scale
		yjitter *= pixel_scale
	
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
		
		if jitter:
			psf = psf.shift(xjitter,yjitter)
		
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
	h['spectrid'] = spectra_id
	
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
		

