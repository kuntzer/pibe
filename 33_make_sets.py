import numpy as np
import glob
from astropy.io import fits

###################################################################################################
# Defining variables

# path to input simulated PSFs 
sim_dir = "output/psf_nonoise"

# Output dir 
out_dir = "output/datasets"

# Number of simulations in the training set
n_img_train = 50#000

# Number of simulations in the calibration/test set
n_img_test = 100#000

###################################################################################################
# Initialisation