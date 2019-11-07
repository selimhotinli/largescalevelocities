#### produces an nbodykit catalog from websky or some other (on dev)
#### contact: selim.hotinli14@imperial.ac.uk
##### Project the field onto a healpix map: 'maps_field'.
from nbodykit.lab import *
from nbodykit import setup_logging
setup_logging() # log output to stdout
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.interpolate import *
import sys, platform, os
from numpy import genfromtxt
import healpy as hp
from IPython.display import display
import six
import warnings
import fitsio
from nbodykit.source.catalog import FITSCatalog
import random
## \user selections
size = 10000000 # size of the catalog
chi_bins = np.linspace(0,5000, 10) # equal-chi bins (where chi is the comoving distance).
# if know, the size of the intendend healpix map that will be used for redshift
# binned velocity or density maps. each halo maps on to a single pixel in total.
# (many halos can map onto the same pixel, of course)
NSIDE = 2**6
## \user selections
NPIX = 12*NSIDE**2
# initialize the map
maps_overdensity = dict()
maps_numbercount = dict()
# Bins to compute transverse perturbations in
chi_bins = np.linspace(0,5000, 10)
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
	maps_numbercount[ii] = np.zeros(NPIX)
# Load the halo cataloge
f=open('/path/to/halo/catalog/e.g./websky/v0.0/halos.pksc')
N=np.fromfile(f,count=3,dtype=np.int32)[0]
catalog=np.fromfile(f,count=N*10,dtype=np.float32)
catalog=np.reshape(catalog,(N,10))
# Select a subset, randomly from all halos
select = random.sample(np.arange(0,len(catalog)), size)
# Cosmological parameters
omegab = 0.049;
omegac = 0.261;
omegam = omegab + omegac
h      = 0.68;
ns     = 0.965;
sigma8 = 0.81
rho_m_0 = 2.775e11*omegam*h**2 # Msun/Mpc^3
c = 3e5
# More cosmology
H0 = 100*h; nz = 100000; z1 = 1e-6
z2 = 5000.0; za = np.linspace(z1,z2,nz)
dz = za[1]-za[0]
H      = lambda z: H0*np.sqrt(omegam*(1+z)**3+1-omegam)
dchidz = lambda z: c/H(z)
chia = np.cumsum(dchidz(za))*dz
# To cover all halos in the catalog...
chia = np.insert(chia,0,0)
za = np.insert(za,0,0)
zofchi = interp1d(chia,za)
# ################## halo locations in Mpc (comoving) coordinates ##################
x = catalog[:,0];
y = catalog[:,1];
z = catalog[:,2]; # Mpc (comoving)
# ################## the comoving radial distance $\chi$ to halos ##################
chi = np.sqrt(x**2 + y**2 + z**2) # Mpc
# ################ the radius of the halos in the catalogue in kpc ################
R  = catalog[:,6] # kpc
# ################## the velocity of the halos in the catalogue ##################
vx = catalog[:,3];
vy = catalog[:,4];
vz = catalog[:,5] # km/sec
# recond (given chi-binning) the corresponding chi-bin of each halo apriori
catind = np.digitize(chi, chi_bins[:-1])
# and the pixel associated with that halo, as explained above
catpix = hp.pixelfunc.vec2pix(NSIDE, x, y, z)
# # convert to mass, comoving distance, radial velocity, redshfit, RA and DEc
M        = 4*np.pi/3.*rho_m_0*R**3    # Msun
# vrad     = (x*vx + y*vy + z*vz) / chi # km/sec
lencat = len(catalog)
cosmo = cosmology.Planck15
# Generate a NBODYKIT catalog from the data
dset = numpy.empty(lencat, dtype=[('Position',
				   ('f8', 3)),
				  ('Mass', 'f8'),
				  ('r', 'f8'),
				  ('t', 'f8'),
				  ('p', 'f8'),
				  ('Velocity', ('f8', 3)),
				  ('vr', 'f8'),
				  ('vp', 'f8'),
				  ('vt', 'f8') ,
				  ('aH', 'f8'),
				  ('pixel_bin','f8'),
				  ('pixel_index','f8')])
r_ = np.sqrt(x**2 + y**2 + z**2)
p_ = np.arctan(y/(x+1e-13))
t_ = np.arccos(z/(r_+1e-13))
dset['r'] = r_;dset['p'] = p_;dset['t'] = t_
dset['Position'] = np.array([x,y,z]).transpose()
dset['Mass'] = M
dset['Velocity'] = np.array([vx,vy,vz]).transpose()
dset['vr'] = np.sin(t_)*np.cos(p_)*vx + np.sin(t_)*np.sin(p_)*vy + np.cos(t_)*vz
dset['vt'] = np.cos(t_)*np.cos(p_)*vx + np.cos(t_)*np.sin(p_)*vy - np.sin(t_)*vz
dset['vp'] = -np.sin(t_)*vx + np.cos(p_)*vy
redshift = zofchi(r_)
f = cosmo.scale_independent_growth_rate(redshift)
aH = (100.*cosmo.efunc(redshift))/(1.+redshift)
dset['aH'] = aH
dset['pixel_bin'] = np.digitize(chi, chi_bins[:-1])
dset['pixel_index'] = hp.pixelfunc.vec2pix(NSIDE, x, y, z)
# write to a FITS file using fitsio so nbodykit can read it
fitsio.write('/write/catalog/here.fits', dset, extname='Data',clobber=True)
