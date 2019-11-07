#### produces an nbodykit catalog from websky or some other (on dev)
#### contact: selim.hotinli14@imperial.ac.uk
##### Project the field onto a healpix map: 'maps_field'.
from nbodykit import *
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
import voperate
##########################
wantx = 1 # what projections of the 3D 3-velocity one should compute?
wanty = 0
wantz = 0
wantr = 1 # radial (line-of-sight)
wantp = 1 # phi (transverse)
wantt = 1 # theta (transverse)
rsd_err = 1 # add redshift space distotion errors?
pze_err = 1 # add photometric redshift errors? (\sigma=(1+z)*0.03)
Nmesh_ = 1000 # Mesh size
BoxSize= 16000 # BoxSize in Mpc (depends on catalog)
NSIDE = 2**7 # healpix map size (should be matched w catalog upon writing)
##########################
cosmo = cosmology.Planck15
# initialize the catalog
cat = FITSCatalog('/rds/general/user/sch14/ephemeral/moving_lens/Websky_Catalog_big.fits', ext='Data')
# Cosmological parameters
omegab = 0.049;
omegac = 0.261;
omegam = omegab + omegac
h      = 0.68;
ns     = 0.965;
sigma8 = 0.81
rho_m_0 = 2.775e11*omegam*h**2 # Msun/Mpc^3
c = 3e5
# Cosmology
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
r = cat['r'].compute()
red = zofchi(r)
# get the angular coordinates of the objects in the catalogue
ra, dec = transform.CartesianToEquatorial(cat['Position'], observer=[0,0,0])
phi   = (np.pi/180.)*ra.compute()
theta = (np.pi/180.)*dec.compute()
st = np.sin(theta)
cp = np.cos(phi)
sp = np.sin(phi)
ct = np.cos(theta)
# spherical unit vectors
r1 = st*cp ## this is theta
r2 = st*sp ## this is phi
r3 = ct    ## this is z
rhat = np.array([r1,r2,r3])
LOS  = rhat.transpose()
rsd = cat['Velocity']/cat['aH'][:,None] * rhat.transpose()
xyz = cat.compute(cat['Position'])
# now re-calculate the positions with the RSD errors.
xyz_rsd = cat['Position'] - rsd
cat['xyz_rsd'] = xyz_rsd
tpr_rsd = transform.CartesianToSky(xyz, cosmo, cat['Velocity'], observer=[0.,0.,0.])
sz = np.random.normal(0, 0.03*(1+red),len(xyz_rsd))
xyz_rsd_pze = transform.SkyToCartesian(tpr_rsd[0], tpr_rsd[1], tpr_rsd[2]*(1.+sz), cosmo, observer=[0.,0.,0.])
cat['xyz_rsd_pze'] = xyz_rsd_pze
# interlaced=True
cat.attrs['Nmesh']=Nmesh_
cat.attrs['BoxSize']=BoxSize
inner_counter = np.array([0])
NPIX = 12*NSIDE**2
maps_overdensity = dict()
maps_numbercount = dict()
# Bins to compute transverse perturbations in
chi_bins = np.linspace(0,5000,20)
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
	maps_numbercount[ii] = np.zeros(NPIX)
zmax = 3000.
mesh_d = cat.to_mesh(compensated=True, window='tsc', position='Position',Nmesh=Nmesh_)
if rsd==1:
	mesh_d_rsd = cat.to_mesh(compensated=True, window='tsc', position='xyz_rsd' ,Nmesh=Nmesh_)
if pze==1:
	mesh_d_rsd_pze = cat.to_mesh(compensated=True, window='tsc', position='xyz_rsd_pze',Nmesh=Nmesh_)
print("starting")
if wantx or wantr or wantp or wantp:
	print("generate x maps")
	maps_overdensity = dict()
	maps_numbercount = dict()
	for ii in range(0,len(chi_bins)):
		maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
	md_obj = mesh_d.apply(get_coefficent,mode='real', kind='relative')
	velx   = md_obj.apply(get_vx, mode='complex', kind='wavenumber')
	velx_  = velx.apply(get_t,mode='real', kind='relative')
	velx__ = velx_.compute()
	maps_velx = maps_overdensity
	maps_velx_counts = maps_numbercount
	for ii in range(0,len(chi_bins)):
		hp.write_map('/path/to/velocity/value/output'+str(ii)+'.fits',maps_velx[ii])
		hp.write_map('/path/to/velocity/mesh/numbercount/output'+str(ii)+'.fits',maps_velx_counts[ii])
	if rsd==1:
		print("generate x maps with rsd")
		maps_overdensity = dict()
		maps_numbercount = dict()
		for ii in range(0,len(chi_bins)):
			maps_overdensity[ii] = np.zeros(NPIX)
			maps_numbercount[ii] = np.zeros(NPIX)
		md_obj = mesh_d_rsd.apply(get_coefficent,mode='real', kind='relative')
		velx   = md_obj.apply(get_vx, mode='complex', kind='wavenumber')
		velx_  = velx.apply(get_t,mode='real', kind='relative')
		velx__ = velx_.compute()
		maps_velx = maps_overdensity
		maps_velx_counts = maps_numbercount
		for ii in range(0,len(chi_bins)):
			hp.write_map('/path/to/velocity/value/output'+str(ii)+'.fits',maps_velx[ii])
			hp.write_map('/path/to/velocity/value/output'+str(ii)+'.fits',maps_velx_counts[ii])
		maps_overdensity = dict()
		maps_numbercount = dict()
	if pze==1:
		print("generate x maps with rsd and pze")
		for ii in range(0,len(chi_bins)):
			maps_overdensity[ii] = np.zeros(NPIX)
			maps_numbercount[ii] = np.zeros(NPIX)
		md_obj = mesh_d_rsd_pze.apply(get_coefficent,mode='real', kind='relative')
		velx   = md_obj.apply(get_vx, mode='complex', kind='wavenumber')
		velx_  = velx.apply(get_t,mode='real', kind='relative')
		velx__ = velx_.compute()
		maps_velx = maps_overdensity
		maps_velx_counts = maps_numbercount
		for ii in range(0,len(chi_bins)):
			hp.write_map('/path/to/velocity/value/output'+str(ii)+'.fits',maps_velx[ii])
			hp.write_map('/path/to/velocity/value/output'+str(ii)+'.fits',maps_velx_counts[ii])
	if wanty or wantr or wantp or wantp:
		print("generate y maps")
		maps_overdensity = dict()
		maps_numbercount = dict()
		for ii in range(0,len(chi_bins)):
			maps_overdensity[ii] = np.zeros(NPIX)
			maps_numbercount[ii] = np.zeros(NPIX)
		md_obj = mesh_d.apply(get_coefficent,mode='real', kind='relative')
		vely   = md_obj.apply(get_vy, mode='complex', kind='wavenumber')
		vely_  = vely.apply(get_t,mode='real', kind='relative')
		vely__ = vely_.compute()
		maps_vely = maps_overdensity
		maps_vely_counts = maps_numbercount
		for ii in range(0,len(chi_bins)):
			hp.write_map('/path/to/velocity/value/output'+str(ii)+'.fits',maps_vely[ii])
			hp.write_map('/path/to/velocity/mesh/numbercount/output'+str(ii)+'.fits',maps_vely_counts[ii])
		if rsd==1:
			print("generate y maps with rsd")
			maps_overdensity = dict()
			maps_numbercount = dict()
			for ii in range(0,len(chi_bins)):
				maps_overdensity[ii] = np.zeros(NPIX)
				maps_numbercount[ii] = np.zeros(NPIX)
			md_obj = mesh_d_rsd.apply(get_coefficent,mode='real', kind='relative')
			vely   = md_obj.apply(get_vy, mode='complex', kind='wavenumber')
			vely_  = vely.apply(get_t,mode='real', kind='relative')
			vely__ = vely_.compute()
			maps_vely = maps_overdensity
			maps_vely_counts = maps_numbercount
			for ii in range(0,len(chi_bins)):
				hp.write_map('/path/to/velocity/value/output'+str(ii)+'.fits',maps_vely[ii])
				hp.write_map('/path/to/velocity/value/output'+str(ii)+'.fits',maps_vely_counts[ii])
			maps_overdensity = dict()
			maps_numbercount = dict()
		if pze==1:
			print("generate y maps with rsd and pze")
			for ii in range(0,len(chi_bins)):
				maps_overdensity[ii] = np.zeros(NPIX)
				maps_numbercount[ii] = np.zeros(NPIX)
			md_obj = mesh_d_rsd_pze.apply(get_coefficent,mode='real', kind='relative')
			vely   = md_obj.apply(get_vy, mode='complex', kind='wavenumber')
			vely_  = vely.apply(get_t,mode='real', kind='relative')
			vely__ = vely_.compute()
			maps_vely = maps_overdensity
			maps_vely_counts = maps_numbercount
			for ii in range(0,len(chi_bins)):
				hp.write_map('/path/to/velocity/value/output'+str(ii)+'.fits',maps_vely[ii])
				hp.write_map('/path/to/velocity/value/output'+str(ii)+'.fits',maps_vely_counts[ii])
	if wantz or wantr or wantp or wantp:
	print("generate z maps")
	maps_overdensity = dict()
	maps_numbercount = dict()
	for ii in range(0,len(chi_bins)):
		maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
	md_obj = mesh_d.apply(get_coefficent,mode='real', kind='relative')
	velz   = md_obj.apply(get_vz, mode='complex', kind='wavenumber')
	velz_  = velz.apply(get_t,mode='real', kind='relative')
	velz__ = velz_.compute()
	maps_velz = maps_overdensity
	maps_velz_counts = maps_numbercount
	for ii in range(0,len(chi_bins)):
		hp.write_map('/path/to/velocity/value/output'+str(ii)+'.fits',maps_velz[ii])
		hp.write_map('/path/to/velocity/mesh/numbercount/output'+str(ii)+'.fits',maps_velz_counts[ii])
	if rsd==1:
		print("generate z maps with rsd")
		maps_overdensity = dict()
		maps_numbercount = dict()
		for ii in range(0,len(chi_bins)):
			maps_overdensity[ii] = np.zeros(NPIX)
			maps_numbercount[ii] = np.zeros(NPIX)
		md_obj = mesh_d_rsd.apply(get_coefficent,mode='real', kind='relative')
		velz   = md_obj.apply(get_vz, mode='complex', kind='wavenumber')
		velz_  = velz.apply(get_t,mode='real', kind='relative')
		velz__ = velz_.compute()
		maps_velz = maps_overdensity
		maps_velz_counts = maps_numbercount
		for ii in range(0,len(chi_bins)):
			hp.write_map('/rds/general/user/sch14/ephemeral/moving_lens/maps_big_1p5KM_AllH_velz_rsd_'+str(ii)+'.fits',maps_velz[ii])
			hp.write_map('/rds/general/user/sch14/ephemeral/moving_lens/maps_big_1p5KM_AllH_velz_rsd_counts_'+str(ii)+'.fits',maps_velz_counts[ii])
		maps_overdensity = dict()
		maps_numbercount = dict()
	if pze==1:
		print("generate z maps with rsd and pze")
		for ii in range(0,len(chi_bins)):
			maps_overdensity[ii] = np.zeros(NPIX)
			maps_numbercount[ii] = np.zeros(NPIX)
		md_obj = mesh_d_rsd_pze.apply(get_coefficent,mode='real', kind='relative')
		velz   = md_obj.apply(get_vz, mode='complex', kind='wavenumber')
		velz_  = velz.apply(get_t,mode='real', kind='relative')
		velz__ = velz_.compute()
		maps_velz = maps_overdensity
		maps_velz_counts = maps_numbercount
		for ii in range(0,len(chi_bins)):
			hp.write_map('/rds/general/user/sch14/ephemeral/moving_lens/maps_big_1p5KM_AllH_velz_'+str(ii)+'.fits',maps_velz[ii])
		hp.write_map('/rds/general/user/sch14/ephemeral/moving_lens/maps_big_1p5KM_AllH_velz_counts_'+str(ii)+'.fits',maps_velz_counts[ii])

maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
velx_theta   = velx.apply(get_v_theta_vx, mode='real', kind='relative')
velx_theta_  = velx_theta.apply(get_t,mode='real', kind='relative')
velx_theta__ = velx_theta_.compute()
maps_velx_theta = maps_overdensity
maps_velx_theta_counts = maps_numbercount
###################
maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
velx_theta_rsd   = velx_rsd.apply(get_v_theta_vx, mode='real', kind='relative')
velx_theta_rsd_  = velx_theta_rsd.apply(get_t,mode='real', kind='relative')
velx_theta_rsd__ = velx_theta_rsd_.compute()
maps_velx_theta_rsd = maps_overdensity
maps_velx_theta_rsd_counts = maps_numbercount
###################
maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
velx_theta_rsd_pze   = velx_rsd_pze.apply(get_v_theta_vx, mode='real', kind='relative')
velx_theta_rsd_pze_  = velx_theta_rsd_pze.apply(get_t,mode='real', kind='relative')
velx_theta_rsd_pze__ = velx_theta_rsd_pze_.compute()
maps_velx_theta_rsd_pze = maps_overdensity
maps_velx_theta_rsd_pze_counts = maps_numbercount
###################
###################
maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
vely_theta   = vely.apply(get_v_theta_vy, mode='real', kind='relative')
vely_theta_  = vely_theta.apply(get_t,mode='real', kind='relative')
vely_theta__ = vely_theta_.compute()
maps_vely_theta = maps_overdensity
maps_vely_theta_counts = maps_numbercount
###################
maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
vely_theta_rsd   = vely_rsd.apply(get_v_theta_vy, mode='real', kind='relative')
vely_theta_rsd_  = vely_theta_rsd.apply(get_t,mode='real', kind='relative')
vely_theta_rsd__ = vely_theta_rsd_.compute()
maps_vely_theta_rsd = maps_overdensity
maps_vely_theta_rsd_counts = maps_numbercount
###################
maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
vely_theta_rsd_pze   = vely_rsd_pze.apply(get_v_theta_vy, mode='real', kind='relative')
vely_theta_rsd_pze_  = vely_theta_rsd_pze.apply(get_t,mode='real', kind='relative')
vely_theta_rsd_pze__ = vely_theta_rsd_pze_.compute()
maps_vely_theta_rsd_pze = maps_overdensity
maps_vely_theta_rsd_pze_counts = maps_numbercount
###################
###################
maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
velz_theta   = velz.apply(get_v_theta_vz, mode='real', kind='relative')
velz_theta_  = velz_theta.apply(get_t,mode='real', kind='relative')
velz_theta__ = velz_theta_.compute()
maps_velz_theta = maps_overdensity
maps_velz_theta_counts = maps_numbercount
###################
maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
velz_theta_rsd   = velz_rsd.apply(get_v_theta_vz, mode='real', kind='relative')
velz_theta_rsd_  = velz_theta_rsd.apply(get_t,mode='real', kind='relative')
velz_theta_rsd__ = velz_theta_rsd_.compute()
maps_velz_theta_rsd = maps_overdensity
maps_velz_theta_rsd_counts = maps_numbercount
###################
maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
velz_theta_rsd_pze   = velz_rsd_pze.apply(get_v_theta_vz, mode='real', kind='relative')
velz_theta_rsd_pze_  = velz_theta_rsd_pze.apply(get_t,mode='real', kind='relative')
velz_theta_rsd_pze__ = velz_theta_rsd_pze_.compute()
maps_velz_theta_rsd_pze = maps_overdensity
maps_velz_theta_rsd_pze_counts = maps_numbercount
###################
maps_theta = dict()
maps_theta_counts = dict()
maps_theta_rsd = dict()
maps_theta_rsd_counts = dict()
maps_theta_rsd_pze = dict()
maps_theta_rsd_pze_counts = dict()
for ii in range(0,10):
	maps_theta[ii] = maps_velx_theta[ii]+maps_vely_theta[ii]+maps_velz_theta[ii]
		maps_theta_counts[ii] = maps_velx_theta_counts[ii]+maps_vely_theta_counts[ii]+maps_velz_theta_counts[ii]
			
			maps_theta_rsd[ii] = maps_velx_theta_rsd[ii]+maps_vely_theta_rsd_pze[ii]+maps_velz_theta_rsd[ii]
    maps_theta_rsd_counts[ii] = maps_velx_theta_rsd_counts[ii]+maps_vely_theta_rsd_counts[ii]+maps_velz_theta_rsd_counts[ii]
		
		maps_theta_rsd_pze[ii] = maps_velx_theta_rsd_pze[ii]+maps_vely_theta_rsd_pze[ii]+maps_velz_theta_rsd_pze[ii]
		maps_theta_rsd_pze_counts[ii] = maps_velx_theta_rsd_pze_counts[ii]+maps_vely_theta_rsd_pze_counts[ii]+maps_velz_theta_rsd_pze_counts[ii]
maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
velx_r   = velx.apply(get_v_r_vx, mode='real', kind='relative')
velx_r_  = velx_r.apply(get_t,mode='real', kind='relative')
velx_r__ = velx_r_.compute()
maps_velx_r = maps_overdensity
maps_velx_r_counts = maps_numbercount
###################
maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
velx_r_rsd   = velx_rsd.apply(get_v_r_vx, mode='real', kind='relative')
velx_r_rsd_  = velx_r_rsd.apply(get_t,mode='real', kind='relative')
velx_r_rsd__ = velx_r_rsd_.compute()
maps_velx_r_rsd = maps_overdensity
maps_velx_r_rsd_counts = maps_numbercount
###################
maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
velx_r_rsd_pze   = velx_rsd_pze.apply(get_v_r_vx, mode='real', kind='relative')
velx_r_rsd_pze_  = velx_r_rsd_pze.apply(get_t,mode='real', kind='relative')
velx_r_rsd_pze__ = velx_r_rsd_pze_.compute()
maps_velx_r_rsd_pze = maps_overdensity
maps_velx_r_rsd_pze_counts = maps_numbercount
###################
###################
maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
vely_r   = vely.apply(get_v_r_vy, mode='real', kind='relative')
vely_r_  = vely_r.apply(get_t,mode='real', kind='relative')
vely_r__ = vely_r_.compute()
maps_vely_r = maps_overdensity
maps_vely_r_counts = maps_numbercount
###################
maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
vely_r_rsd   = vely_rsd.apply(get_v_r_vy, mode='real', kind='relative')
vely_r_rsd_  = vely_r_rsd.apply(get_t,mode='real', kind='relative')
vely_r_rsd__ = vely_r_rsd_.compute()
maps_vely_r_rsd = maps_overdensity
maps_vely_r_rsd_counts = maps_numbercount
###################
maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
vely_r_rsd_pze   = vely_rsd_pze.apply(get_v_r_vy, mode='real', kind='relative')
vely_r_rsd_pze_  = vely_r_rsd_pze.apply(get_t,mode='real', kind='relative')
vely_r_rsd_pze__ = vely_r_rsd_pze_.compute()
maps_vely_r_rsd_pze = maps_overdensity
maps_vely_r_rsd_pze_counts = maps_numbercount
###################
###################
maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
velz_r   = velz.apply(get_v_r_vz, mode='real', kind='relative')
velz_r_  = velz_r.apply(get_t,mode='real', kind='relative')
velz_r__ = velz_r_.compute()
maps_velz_r = maps_overdensity
maps_velz_r_counts = maps_numbercount
###################
maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
velz_r_rsd   = velz_rsd.apply(get_v_r_vz, mode='real', kind='relative')
velz_r_rsd_  = velz_r_rsd.apply(get_t,mode='real', kind='relative')
velz_r_rsd__ = velz_r_rsd_.compute()
maps_velz_r_rsd = maps_overdensity
maps_velz_r_rsd_counts = maps_numbercount
###################
maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
velz_r_rsd_pze   = velz_rsd_pze.apply(get_v_r_vz, mode='real', kind='relative')
velz_r_rsd_pze_  = velz_r_rsd_pze.apply(get_t,mode='real', kind='relative')
velz_r_rsd_pze__ = velz_r_rsd_pze_.compute()
maps_velz_r_rsd_pze = maps_overdensity
maps_velz_r_rsd_pze_counts = maps_numbercount
###################
maps_r = dict()
maps_r_counts = dict()
maps_r_rsd = dict()
maps_r_rsd_counts = dict()
maps_r_rsd_pze = dict()
maps_r_rsd_pze_counts = dict()

for ii in range(0,10):
	maps_r[ii] = maps_velx_r[ii]+maps_vely_r[ii]+maps_velz_r[ii]
		maps_r_counts[ii] = maps_velx_r_counts[ii]+maps_vely_r_counts[ii]+maps_velz_r_counts[ii]
			
			maps_r_rsd[ii] = maps_velx_r_rsd[ii]+maps_vely_r_rsd_pze[ii]+maps_velz_r_rsd[ii]
    maps_r_rsd_counts[ii] = maps_velx_r_rsd_counts[ii]+maps_vely_r_rsd_counts[ii]+maps_velz_r_rsd_counts[ii]
		
		maps_r_rsd_pze[ii] = maps_velx_r_rsd_pze[ii]+maps_vely_r_rsd_pze[ii]+maps_velz_r_rsd_pze[ii]
		maps_r_rsd_pze_counts[ii] = maps_velx_r_rsd_pze_counts[ii]+maps_vely_r_rsd_pze_counts[ii]+maps_velz_r_rsd_pze_counts[ii]

maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
velx_phi   = velx.apply(get_v_phi_vx, mode='real', kind='relative')
velx_phi_  = velx_phi.apply(get_t,mode='real', kind='relative')
velx_phi__ = velx_phi_.compute()
maps_velx_phi = maps_overdensity
maps_velx_phi_counts = maps_numbercount
###################
maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
velx_phi_rsd   = velx_rsd.apply(get_v_phi_vx, mode='real', kind='relative')
velx_phi_rsd_  = velx_phi_rsd.apply(get_t,mode='real', kind='relative')
velx_phi_rsd__ = velx_phi_rsd_.compute()
maps_velx_phi_rsd = maps_overdensity
maps_velx_phi_rsd_counts = maps_numbercount
###################
maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
velx_phi_rsd_pze   = velx_rsd_pze.apply(get_v_phi_vx, mode='real', kind='relative')
velx_phi_rsd_pze_  = velx_phi_rsd_pze.apply(get_t,mode='real', kind='relative')
velx_phi_rsd_pze__ = velx_phi_rsd_pze_.compute()
maps_velx_phi_rsd_pze = maps_overdensity
maps_velx_phi_rsd_pze_counts = maps_numbercount
###################
###################
maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
vely_phi   = vely.apply(get_v_phi_vy, mode='real', kind='relative')
vely_phi_  = vely_phi.apply(get_t,mode='real', kind='relative')
vely_phi__ = vely_phi_.compute()
maps_vely_phi = maps_overdensity
maps_vely_phi_counts = maps_numbercount
###################
maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
vely_phi_rsd   = vely_rsd.apply(get_v_phi_vy, mode='real', kind='relative')
vely_phi_rsd_  = vely_phi_rsd.apply(get_t,mode='real', kind='relative')
vely_phi_rsd__ = vely_phi_rsd_.compute()
maps_vely_phi_rsd = maps_overdensity
maps_vely_phi_rsd_counts = maps_numbercount
###################
maps_overdensity = dict()
maps_numbercount = dict()
for ii in range(0,len(chi_bins)):
	maps_overdensity[ii] = np.zeros(NPIX)
		maps_numbercount[ii] = np.zeros(NPIX)
vely_phi_rsd_pze   = vely_rsd_pze.apply(get_v_phi_vy, mode='real', kind='relative')
vely_phi_rsd_pze_  = vely_phi_rsd_pze.apply(get_t,mode='real', kind='relative')
vely_phi_rsd_pze__ = vely_phi_rsd_pze_.compute()
maps_vely_phi_rsd_pze = maps_overdensity
maps_vely_phi_rsd_pze_counts = maps_numbercount
###################
maps_phi = dict()
maps_phi_counts = dict()
maps_phi_rsd = dict()
maps_phi_rsd_counts = dict()
maps_phi_rsd_pze = dict()
maps_phi_rsd_pze_counts = dict()

for ii in range(0,10):
	maps_phi[ii] = maps_velx_phi[ii]+maps_vely_phi[ii]
		maps_phi_counts[ii] = maps_velx_phi_counts[ii]+maps_vely_phi_counts[ii]
			
			maps_phi_phisd[ii] = maps_velx_phi_rsd[ii]+maps_vely_phi_rsd_pze[ii]
    maps_phi_rsd_counts[ii] = maps_velx_phi_rsd_counts[ii]+maps_vely_phi_rsd_counts[ii]
		
		maps_phi_rsd_pze[ii] = maps_velx_phi_rsd_pze[ii]+maps_vely_phi_rsd_pze[ii]
		maps_phi_rsd_pze_counts[ii] = maps_velx_phi_rsd_pze_counts[ii]+maps_vely_phi_rsd_pze_counts[ii]

