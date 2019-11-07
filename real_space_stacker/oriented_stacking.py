#### produces code to be updated and under development
#### contact: selim.hotinli14@imperial.ac.uk
#### produce stacked imagines of CMB around halos
import sys, platform, os
from sys import argv
import numpy as np
from numpy import genfromtxt
from matplotlib import pyplot as plt
import healpy as hp
import scipy
from scipy.interpolate import *
from IPython.display import display
import six
import warnings
c = 3e5
# Set your cosmology below
omegab = 0.049
omegac = 0.261
omegam = omegab + omegac
h      = 0.68
ns     = 0.965
sigma8 = 0.81
H0 = 100*h
nz = 100000
z1 = 0.0
z2 = 100.0
za = np.linspace(z1,z2,nz)
dz = za[1]-za[0]
H      = lambda z: H0*np.sqrt(omegam*(1+z)**3+1-omegam)
dchidz = lambda z: c/H(z)
chia = np.cumsum(dchidz(za))*dz
zofchi = interp1d(chia,za)
# Load in catalog
f=open('/path/to/catalog/e.g./websky/v0.0/halos.pksc')
N=np.fromfile(f,count=3,dtype=np.int32)[0]
catalog=np.fromfile(f,count=N*10,dtype=np.float32)
catalog=np.reshape(catalog,(N,10))
# Mpc (comoving)
x = catalog[:,0];
y = catalog[:,1];
z = catalog[:,2];
chi_ = np.sqrt(x**2 + y**2 + z**2) # Mpc
# Can bin the data for memory saving and/or speed when run in parallel
chi_bins = np.linspace(0,5000,200)
for bin_num in range(len(chi_bins)-1) :	
	if bin_num != int(argv[1]): continue
	print("Working in bin", bin_num, "of", len(chi_bins), "...")
	print("Reducing catalogue...")
    	bin_lower_chi = chi_bins[bin_num]
	bin_upper_chi = chi_bins[bin_num+1]
	d_chi = bin_upper_chi - bin_lower_chi
	reduced_catalog = catalog[ (bin_lower_chi<chi_) & (chi_<bin_upper_chi) ]
	print("...done.")
	cutCat = np.arange(0,len(catalog))
	cutCat = cutCat[(bin_lower_chi<chi_) & (chi_<bin_upper_chi)]		
	NSIDE = 4096
	NPIX = 12*NSIDE**2
	rho_m_0 = 2.775e11*omegam*h**2 # Msun/Mpc^3
	################## halo locations in Mpc (comoving) coordinates ##################
	x = reduced_catalog[:,0]
	y = reduced_catalog[:,1]
	z = reduced_catalog[:,2] # Mpc (comoving)
	################## the comoving radial distance $\chi$ to halos ##################
	chi = np.sqrt(x**2 + y**2 + z**2) # Mpc
	################ the radius of the halos in the catalogue in kpc ################
	R  = reduced_catalog[:,6] # kpc
	################## the velocity of the halos in the catalogue ###################
	vx = reduced_catalog[:,3]
	vy = reduced_catalog[:,4]
	vz = reduced_catalog[:,5] # km/sec
	pix = hp.vec2pix(NSIDE, x, y, z)
	### convert to mass, comoving distance, radial velocity, redshfit, RA and DEc ###
	M        = 4*np.pi/3.*rho_m_0*R**3    # Msun
	vrad     = (x*vx + y*vy + z*vz) / chi # km/sec
	redshift = zofchi(chi)
	# Below is a useful quantity: scale radious 'srad' as defined in the literature.
	## halo-model dependent constants, can be found in literature ##
	Ac = 7.85
	alphac = -0.081
	betac  = -0.71
	cs = Ac*(M/(2e12/h))**alphac*(1.+redshift)**betac
	srad = R/cs  # Mpc
	## The angular values corresponding to the sky locations of the halos ##
	theta, phi  = hp.vec2ang(np.column_stack((x,y,z))) # in radians
	## The dimensionless number useful in identifying the size of the dipolar modulation ##
	xfactors = chi/srad /(1+redshift)
	############################################################
	# This is simply the orientation of the center of the patch
	#  that's already been carried on to the middle of the map.
	mid_a = hp.ang2pix(NSIDE,np.pi/2.,0)
	mid_v = hp.pix2vec(NSIDE,mid_a)
	x0 = mid_v[0];
	y0 = mid_v[1];
	z0 = mid_v[2];
	############################################################
		

	############################################################
	# ___ We calculate the angular distances that amounts to a certain distance in radii.
	# 1 ) choose a window size as some number times the angular size corresponding to the
	#     scale-radious at the objects distance (in redshift).
	# 1.! ) make flat sky approximation, i.e. radious/distance = angle.
	window = 6.
	#
	# srad is the object's scale radious and chi is the distance from the observer, in comoving units.
	srad_over_chi = srad/chi
	angular_dist = srad_over_chi*window # angular diameter distance
	#####################################################
	#### Calculate the transfer matrix for all halos ####
	#####################################################
	xhat = x/chi;
	yhat = y/chi;
	zhat = z/chi # unit distance vectors
	v_rad = np.sqrt(vx**2.+vy**2.+vz**2.) # velocity amplitude of all halos (km / s)
	vtheta = ( np.cos(theta)*np.cos(phi)*vx + np.cos(theta)*np.sin(phi)*vy - np.sin(theta)*vz )
	vphi   = (-np.sin(phi)*vx + np.cos(phi)*vy )
	thetahat = np.array([np.cos(theta)*np.cos(phi),np.cos(theta)*np.sin(phi), - np.sin(theta)])
	phihat   = np.array([ -np.sin(phi),np.cos(phi),0.*phi])
	# transverse velocity vector
	vtrans_ = vtheta*thetahat + vphi*phihat
	vtans_x = vtrans_[0,:]
	vtans_y = vtrans_[1,:]
	vtans_z = vtrans_[2,:]
	vtans_r = np.sqrt(vtans_x**2.+vtans_y**2.+vtans_z**2.)
	############# the transformation matrix #############
	#
	#      i_1 = r
	#      i_2 = v_t
	#      i_3 = r x v_t
	#
	#        A_1'    i_1'*i_1  i_1'*i_2  i_1'*i_3    A_1
	#        A_2'  = i_2'*i_1  i_2'*i_2  i_2'*i_3 *  A_2
	#        A_3'    i_3'*i_1  i_3'*i_2  i_3'*i_3    A_3
	#      _output_  ____________\/______________  _input_
	#                          mTR
	#      i_1' = x
	#      i_2' = y
	#     i_3' = z
	#
	####################################################
	i1_v = np.transpose([xhat,yhat,zhat])
	i2_x = vtans_x/vtans_r
	i2_y = vtans_y/vtans_r
	i2_z = vtans_z/vtans_r
	i2_v = np.transpose([i2_x,i2_y,i2_z])
		
	i3_x = yhat * i2_z - zhat * i2_y
	i3_y = zhat * i2_x - xhat * i2_z
	i3_z = xhat * i2_y - yhat * i2_x
	i3_v = np.transpose([i3_x,i3_y,i3_z])
		
	i1_v_p = np.array([1.,0.,0.])
	i2_v_p = np.array([0.,1.,0.])
	i3_v_p = np.array([0.,0.,1.])
		
	mTR = np.empty((len(i3_v),3,3),dtype=float)
		
	for ii in range(0,len(i3_v)):
		mTR[ii] = np.array([[np.dot(i1_v_p,i1_v[ii]),np.dot(i1_v_p,i2_v[ii]),np.dot(i1_v_p,i3_v[ii])],
												[np.dot(i2_v_p,i1_v[ii]),np.dot(i2_v_p,i2_v[ii]),np.dot(i2_v_p,i3_v[ii])],
												[np.dot(i3_v_p,i1_v[ii]),np.dot(i3_v_p,i2_v[ii]),np.dot(i3_v_p,i3_v[ii])]])
	# if the radii are given as comoving distance then this amounts to comoving angular distance?
	theta_comov = srad_over_chi*4.
	ang_res_NSIDE = hp.pixelfunc.nside2resol(NSIDE) #
	SO_res_1p4amn = 1.4/3437.75
	npixelsSO = theta_comov/SO_res_1p4amn
	npixels4K = theta_comov/ang_res_NSIDE
	cmb_map  = hp.read_map('/path/to/CMB/map')
	movl_map = hp.read_map('/path/to/moving-lens/map')
	print("read the maps")
	cmb_cls  = np.loadtxt('/path/to/CMB/spectra')
	movl_cls = np.loadtxt('/path/to/moving-lens/spectra')
	print("read the cls")
	### correlation functions if needed (currently disabled) ###
	#			cth_cmb = np.loadtxt('/path/to/CMB/correlation-function')
	#			cth_mov = np.loadtxt('/path/to/moving-lens/correlation-function')
	#			cth_tot = np.loadtxt('/path/to/total CMB (inc. moving lens)/correlation-function')
	### correlation functions if needed (currently disabled) ###
	# 		mu = np.linspace(0.5,1,100000)
	#     th = np.arccos(mu)
	#     corr_cmb = interpolate.interp1d(th, cth_cmb)
	#     corr_mov = interpolate.interp1d(th, cth_mov)
	#     corr_tot = interpolate.interp1d(th, cth_tot)
	print("read the corr funcs")
	ell = np.linspace(2,10000,9999)
	fac1 = 1/(8.*np.log(2.))
	fac2 = ell * (ell + 1. ) * fac1
	# various experimental RMS and beam details
	fwhmSO = 1.4/3437.75
	deltSO = 6./3437.75
	fwhmS4 = 1.4/3437.75
	deltS4 = 1./3437.75
	facs = fwhmSO**2
	facd = deltSO**2
	facsS4 = fwhmS4**2
	facdS4 = deltS4**2
	TCMB = 2.725 # in Kelvin
	muK2 = (TCMB*1e6)**2
	# idealized CMB noise
	NTT   = facd*np.exp(fac2 * facs) #* (TCMB*1e6)**2 # in [mK]^2
	NTTS4 = facdS4*np.exp(fac2 * facsS4) #* (TCMB*1e6)**2 # in [mK]^2
	
	# Now, we wish to create a template on-to which we will *map* our CMB patch.
	# This is a template in dimensionless units 'x' as described in the notes,
	# as well as above. We should decide what range of 'x' we should consider
	# as the signal profile is largely captured by around x=2.
		
	# Choice of grid-size and binning is arbitrary
	x_h = np.linspace(-3., 3., 11,dtype=float)
	y_h = np.linspace(-3., 3., 11,dtype=float)
	xv, yv = np.meshgrid(x_h, y_h)		
	xvv = np.reshape(xv, len(xv)*len(xv[0]))
	yvv = np.reshape(yv, len(yv)*len(yv[0]))		
	indvv = range(0,len(xvv))
	leni = len(indvv)		
	numindices = np.array([])
	cmb_map  = unl_cmb_map_4K_f1p4
	mov_map  = movl_map*np.sqrt(muK2)
	tot_map  = mov_map + cmb_map
	stack_s = len(catalog)
	##########
	grid_data  = np.zeros((3,leni),dtype=float)
	##########
	lentmp = len(xvv)
	counter = 0
	resol = hp.nside2resol(NSIDE)
	iselect = np.array([])
	datacmb = open('/path/to/output/grid/per/redshift/bin/'+('%03d' % int(argv[1]))+'.txt', 'w')
	datamov = open('/path/to/output/grid/per/redshift/bin/'+('%03d' % int(argv[1]))+'.txt', 'w')
	datatot = open('/path/to/output/grid/per/redshift/bin/'+('%03d' % int(argv[1]))+'.txt', 'w')
	for index, ii in enumerate(cutCat):
	########################################
	# OK SO HERE: We are rescaling our *grid* to fit well onto
	# ----------  the (expected) shape of the moving-lens effect.
	# ----------  Below are theta and phi-directed template coordinates.
	########################################
		x_t =  yvv/xfactors[index]+np.pi/2.
		x_p = -xvv/xfactors[index]
		# get the (x,y,z) vector corresponding to the *template* coordinates
		vec = hp.ang2vec(x_t,x_p)
		x_ = vec[:,0]; y_ = vec[:,1]; z_ = vec[:,2]
		iselect = np.append(iselect,index)
		patch = np.empty((len(y_),3),dtype=float)
		for jj in range(0,len(x_)):
			patch[jj] = hp.rotator.rotateVector(mTR[index],x_[jj],vy=y_[jj],vz=z_[jj],do_rot=True)
			newp = hp.vec2pix(NSIDE,patch[:,0],patch[:,1],patch[:,2])
			grid_data[0] += tot_map[newp] - tot_map[60]
			grid_data[1] += cmb_map[newp] - cmb_map[60]
			grid_data[2] += mov_map[newp] - mov_map[60]
			if(np.mod(counter,500)==0): print("counter = ", counter)
			counter = counter + 1
		np.savetxt(datacmb, grid_data[1] , delimiter=',', newline="\n")
		np.savetxt(datamov, grid_data[2] , delimiter=',', newline="\n")
		np.savetxt(datatot, grid_data[0] , delimiter=',', newline="\n")
	datacmb.close()
	datamov.close()
	datatot.close()
