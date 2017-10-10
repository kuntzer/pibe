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
sed = 3

cat_dir = "output/psf_nonoise_smallpx_SED{}".format(sed)

out_dir = "output/interpolations"
out_name = "interp_SED{}.pkl".format(sed)

plot_white = False
show = False

###################################################################################################
# Initialisation
cat = Table.read(os.path.join(cat_dir, "catsimmeas.fits"))

print "Using {} simulated PSF shape samples".format(len(cat))

###################################################################################################
# Interpolating the PSF field

x = np.array(cat['x'])
y = np.array(cat['y'])
interp_g1 = interp.SmoothBivariateSpline(x, y, cat['g1'])
interp_g2 = interp.SmoothBivariateSpline(x, y, cat['g2'])
interp_fwhm = interp.SmoothBivariateSpline(x, y, cat['fwhm'])


a = np.linspace(0, 1, 25)
b = np.linspace(0, 1, 25)
c, d = np.meshgrid(a, b)

c = c.flatten()
d = d.flatten()

if plot_white:
	import matplotlib.font_manager as fm
	from matplotlib import rc
	import matplotlib
	prop = fm.FontProperties(fname='/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyreadventor-regular.otf')
	rc('font', **{'family':'TeX Gyre Adventor','size':14})
	matplotlib.rcParams.update({'text.color': 'white', 'ytick.color':'white', 'xtick.color':'white',
		'axes.labelcolor':'white'})
	
	fig = plt.figure()
	ax = plt.subplot(111)
	ax.spines['bottom'].set_color('white')
	ax.spines['top'].set_color('white') 
	ax.spines['right'].set_color('white')
	ax.spines['left'].set_color('white')
	plt.quiver(c, d, interp_g1.ev(c, d), interp_g2.ev(c, d), interp_fwhm.ev(c, d)/10., cmap="Oranges")
	plt.xlabel("x field [deg]")
	plt.ylabel("y field [deg]")
	cb = plt.colorbar()
	cb.set_label("FWHM ['']")
	
	fig.savefig("pibe_psfs.pdf", transparent=True)
	plt.show()
	exit()

# Plotting for visual inspection
for intp in ['g1', 'g2', 'fwhm']: 
	plt.figure()
	a = np.linspace(0, 1, 100)
	b = np.linspace(0, 1, 100)
	c, d = np.meshgrid(a, b)
	plt.title(intp)
	cb = plt.imshow(eval("interp_{}.ev(c,d)".format(intp)), interpolation="None")
	plt.colorbar(cb)
	plt.xlabel("field of view x-axis")
	plt.ylabel("field of view y-axis")
	if show:
		plt.show()
	
# Now saving!
intfn = os.path.join(out_dir, out_name)
utils.writepickle([interp_g1, interp_g2, interp_fwhm], intfn)
print 'Saved interpolants to {}'.format(intfn)
	
