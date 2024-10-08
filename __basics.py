import matplotlib

matplotlib.use('agg')
from scipy.interpolate import interp1d, interp2d, interpn, CubicSpline
from skbio.stats.composition import clr_inv
from astropy import constants as const
import matplotlib.pyplot as plt
from spectres import spectres
import scipy as sp
import numpy as np
import platform
import scipy.io
import warnings
import random
import json
import time
import glob
import math
import copy
import sys
import os

np.set_printoptions(threshold=2**31-1)
warnings.filterwarnings('ignore')
np.random.seed(42)  # DON'T PANIC !!
