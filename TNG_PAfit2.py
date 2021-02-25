#Estimate the kinematic position angle for gas, stars and dark matter at 2Re Illustris-TNG simualtion
#----------------------------------------------------------------------------------------------------
#--------------------------------------------mojtabaraouf@gmail.com----------------------------------
#TNG project-----------------------------------------------------------------------------------------
## Release soon(do not hesitate to send me an email----------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from scipy import stats
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
from   pafit.fit_kinematic_pa import fit_kinematic_pa
# start up your interface of choice and define a helper function, whose purpose is to make a HTTP GET request to a specified URL ("endpoint"), and verify that the response is successful.
#---------------------------------------------------------------------------------------------
import requests
