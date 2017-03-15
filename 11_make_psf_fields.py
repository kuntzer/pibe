from astropy.table import Table
from matplotlib.patches import Ellipse
import pylab as plt
import scipy.interpolate as interp
import os
import numpy as np 

import utils 

###################################################################################################
# Defining variables

# Directory containing psf measurments
cat_dir = "inputdata/fullpsffield"

###################################################################################################
# Initialisation
cat = Table.read(os.path.join(cat_dir, "catmeas.fits"))

print "Using {} PSF shape samples".format(len(cat))

utils.colnorm(cat, "xfield", "x")
utils.colnorm(cat, "yfield", "y")

###################################################################################################
# Interpolating the PSF field

x = np.array(cat['x'])
y = np.array(cat['y'])
interp_g1 = interp.SmoothBivariateSpline(x, y, cat['g1'])
interp_g2 = interp.SmoothBivariateSpline(x, y, cat['g2'])
interp_fwhm = interp.SmoothBivariateSpline(x, y, cat['fwhm'])

# Plotting for visual inspection
for intp in ['g1', 'g2', 'fwhm']: 
	plt.figure()
	a = np.linspace(0, 1, 100)
	b = np.linspace(0, 1, 100)
	c, d = np.meshgrid(a, b)
	plt.title(intp)
	cb = plt.imshow(eval("interp_{}.ev(c,d)".format(intp)), interpolation="None")
	plt.colorbar(cb)
	plt.show()
	
# Now saving!
intfn = os.path.join(cat_dir, "interp.pkl")
utils.writepickle([interp_g1, interp_g2, interp_fwhm], intfn)
print 'Saved interpolants to {}'.format(intfn)
	
