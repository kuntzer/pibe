import os
import glob
import numpy as np
from astropy.io import fits
from astropy.table import Table
import ntpath
from shutil import copyfile


import pylab as plt

###################################################################################################
# Defining variables

# Dataset dir 
dset_dir = "output/datasets_bigpx"

# Output dir 
out_dir = "output/datasets_packaged"

target_snr = [20, 100, -1] 

debug_show = False

###################################################################################################
# initialisation

names = ["train", "test", "calib"]

if not os.path.exists(out_dir):
	os.mkdir(out_dir)
	
target_snr = np.array(target_snr)
sigmas = 1. / (np.pi * (1.25 * 1.225)**2 * target_snr)

###################################################################################################
# Engine start 

for name in names:
	
	print "Starting on data set {}".format(name)
	
	for sigma, tsnr in zip(sigmas, target_snr):
		sigma_ori = sigma
		
		print "\tTarget SNR is {}, Sigma is {}".format(tsnr, sigma)
		
		iin = os.path.join(dset_dir, name)
		cout = os.path.join(out_dir, name)
		
		if not os.path.exists(cout):
			os.mkdir(cout)
			
		catdirout = os.path.join(out_dir, name, "catalogs")
		
		if not os.path.exists(catdirout):
			os.mkdir(catdirout)
		
		if tsnr < 0:	
			cout = os.path.join(out_dir, name, "snr_no")
		else:
			cout = os.path.join(out_dir, name, "snr_{}".format(int(tsnr)))
		
		if not os.path.exists(cout):
			os.mkdir(cout)
			
		imfnames = glob.glob(os.path.join(iin, "{}*[0-9].fits".format(name)))
		
		for ifn, fn in enumerate(imfnames):
			print "\t\tReading {} (img {}/{})".format(fn, ifn+1, len(imfnames))
			ii_fn = int((fn.split(".fits")[0]).split("/{}_".format(name))[1])			
			im = fits.getdata(fn)
			header = fits.getheader(fn)
			
			if tsnr > 0:
				header["tarSNR"] = tsnr
				header["sigma"] = sigma_ori
					
				sigma = sigma_ori * np.ones((header["NLINES"], header["NCOLS"]))
				sigma *=  np.random.normal(1., 0.1, size=sigma.shape)
				
				sigma_table = sigma
				
				sigma = np.repeat(sigma, header["NXPSF1"], 1)
				sigma = np.repeat(sigma, header["NXPSF2"], 0)
				
				cat = Table.read(os.path.join(iin, "{}_{:03}_truth_cat.fits".format(name,ii_fn)))
				
				nstar = cat["fwhm"].flatten().size
				sigma_table = sigma_table.flatten()[:nstar]
				
				snr = 1. / (sigma_table * (1.22*cat["fwhm"])**2 * np.pi)#.reshape((header["NLINES"], header["NCOLS"])))
				snr = snr.flatten() 
	
				im += np.random.normal(scale=sigma, size=im.shape)
				
				print "\t\tMean SNR = {:.2f} +/- {:.2f}".format(snr.mean(), snr.std())
			
			fits.writeto(os.path.join(cout, "train_{:03}.fits".format(ii_fn)), im, header, clobber=True)
			
			os.path.join(out_dir, name)
			fns = glob.glob(os.path.join(iin, "{}_{:03}_*_cat.fits".format(name,ii_fn)))
			for fn in fns:
				copyfile(fn, os.path.join(out_dir, name, "catalogs", ntpath.basename(fn)))
				
			if debug_show:
				plt.figure()
				plt.imshow(im, interpolation="None", cmap=plt.get_cmap("gray"))
				
				if tsnr > 0:
					plt.figure()
					plt.hist(snr,50)
				plt.show()
		
	
	