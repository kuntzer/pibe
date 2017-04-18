import numpy as np
import glob
from astropy.io import fits
from astropy.table import Table
import os

import utils

###################################################################################################
# Defining variables

# path to input simulated PSFs 
sim_dir = "output/psf_nonoise_colour*"

# Interpolation dir 
interp_dir = "output/psf_nonoise_smallpx"

# Output dir 
out_dir = "output/datasets_bigpx_colour"

# Number of simulations in the training set
n_img_train = 196#000

# Number of simulations in the test set
n_img_test = 4#000

# Number of simualtions in the calibration set
n_img_calib = 0#

# Number of maximum psf images
n_max_per_img = 196

###################################################################################################
# Initialisation

iitrain = int(np.ceil(np.sqrt(n_max_per_img)))
iitest = int(np.ceil(np.sqrt(n_max_per_img)))
iicalib = int(np.ceil(np.sqrt(n_max_per_img)))

all_psfs = glob.glob(os.path.join(sim_dir, "star_*.fits"))

set_names = ['train', 'test', 'calib']

if not os.path.exists(os.path.join(out_dir)):
	os.mkdir(os.path.join(out_dir))
	
intfn = os.path.join(interp_dir, "interp.pkl")
interp_g1, interp_g2, interp_fwhm = utils.readpickle(intfn)

used_fns = []
count_img = 0
		
###################################################################################################
# Creating the images + truth & source catalogues

for sn in set_names:
	
	print "Starting on data set {}".format(sn)
		
	icurrent = 0
	ntot = eval("n_img_{}".format(sn))
	iis = eval("ii{}".format(sn))
	
	apfss = all_psfs[icurrent:icurrent+ntot]
	nimgs = int(np.ceil(len(apfss) * 1.0 / n_max_per_img))
	
	if not os.path.exists(os.path.join(out_dir, sn)):
		os.mkdir(os.path.join(out_dir, sn))

	for iimg in range(nimgs):	
		
		print "\tImage {}/{}, data set {}".format(iimg+1, nimgs, sn)
		
		ix = 0
		iy = 0
		
		naxis1 = None
		naxis2 = None
		npxim1 = None
		npxim2 = None
		
		xs = []
		ys = []
		ixs = []
		iys = []
		colours = []
		
		imfn = os.path.join(out_dir, sn, "{}_{:03d}.fits".format(sn, iimg))
		selected_psfs = apfss[iimg * n_max_per_img: (iimg + 1) * n_max_per_img]
		for kkimg, fn_psf in enumerate(selected_psfs):
			used_fns.append(int((fn_psf.split('.fits')[0]).split("/star_")[1]))
			count_img += 1
			
			if kkimg % 10 == 0:
				print "\t\tLoading psf {}/{}".format(kkimg + 1, len(selected_psfs))
				
			psf = fits.getdata(fn_psf)
			flux = psf.sum()
			psf /= flux
			header = fits.getheader(fn_psf)
			
			xs.append(header['X'])
			ys.append(header['Y'])
			colours.append(header['spectrid'])
			
			if naxis1 is None or naxis2 is None:
				naxis1 = header['NAXIS2']
				naxis2 = header['NAXIS1']
				npxim1 = naxis1 * iis
				npxim2 = naxis2 * int(np.ceil((len(selected_psfs) + 0.0) / iis))
				setimg = np.zeros([npxim2, npxim1])
			
			setimg[iy*naxis2:(iy+1)*naxis2, ix * naxis1:(ix+1)*naxis1] = psf

			ixs.append(ix)
			iys.append(iy)
			
			ix += 1
			if ix >= iis:
				ix = 0
				iy += 1
				
		h = fits.Header()
		h["nlines"] = int(np.ceil((len(selected_psfs) + 0.0) / iis))
		h["ncols"] = iis
		h["nxpsf1"] = naxis1
		h["nxpsf2"] = naxis2
		fits.writeto(imfn, setimg, clobber=True, header=h)
		xs = np.array(xs)
		ys = np.array(ys)
		ixs = np.array(ixs)
		iys = np.array(iys)
		
		g1s = interp_g1.ev(xs, ys)
		g2s = interp_g2.ev(xs, ys)
		fwhms = interp_fwhm.ev(xs, ys)
		
		catt = Table([ixs.tolist(), iys.tolist(), ((ixs + 0.5) * naxis1).tolist(), ((iys + 0.5) * naxis2).tolist(),
					((iys + 0.5) * naxis2).tolist(),  ((ixs + 0.5) * naxis1).tolist(), 
					xs, ys, g1s, g2s, fwhms, colours], names=('xcol', 'ycol', 'xcat', 'ycat', 'xpycat', 'ypycat', 'xfield', 'yfield', 'g1', 'g2', 'fwhm', "spectrid"))
		catt['fwhm'].unit = 'arcsec'
		
		cats = Table([ixs.tolist(), iys.tolist(), ((ixs + 0.5) * naxis1).tolist(), ((iys + 0.5) * naxis2).tolist(), 
					((iys + 0.5) * naxis2).tolist(), ((ixs + 0.5) * naxis1).tolist(), 
					xs, ys], 
					names=('xcol', 'ycol', 'xcat', 'ycat', 'xpycat', 'ypycat', 'xfield', 'yfield'))
		
		catfnt = os.path.join(out_dir, sn, "{}_{:03d}_truth_cat.fits".format(sn, iimg))
		catfns = os.path.join(out_dir, sn, "{}_{:03d}_source_cat.fits".format(sn, iimg))
		
		catt.write(catfnt, overwrite=True)
		cats.write(catfns, overwrite=True)
		
		assert np.unique(len(used_fns)) == count_img

	icurrent += ntot

