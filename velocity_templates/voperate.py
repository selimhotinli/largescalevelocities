#### NBODYKIT functions that are used to produce 3D fields.
#### contact: selim.hotinli14@imperial.ac.uk
##### Project the field onto a healpix map: 'maps_field'.
##### -- Counts the total number of mesh points that
#####    contributed to that pixel.
##### -- TODO: Do the integral over the redshift bin better
def get_t(x, v):
	x0 = x[0] + 1e-11
	y0 = x[1]
	z0 = x[2]
	r = np.sqrt(x0**2.+y0**2.+z0**2.)
	# debugging
	x_ = x0*1.+y0*0.+z0*0.
	y_ = x0*0.+y0*1.+z0*0.
	z_ = x0*0.+y0*0.+z0*1.
	p = np.arctan(y_/x_)
	t = np.arccos(z_/r)
	tot_l = len(r)*len(r)
	r = r.reshape(tot_l)
	p = p.reshape(tot_l)
	t = t.reshape(tot_l)
	if x0[0,0] > 0:
		p = np.pi - p
	pixels = hp.pixelfunc.vec2pix(NSIDE, x_, y_, z_)
	pixels = pixels.reshape(tot_l)
	values = v.reshape(tot_l)
	inds = np.digitize(r, chi_bins[:-1])
	for ii in range(0,len(pixels)):
		maps_field[inds[ii]-1][pixels[ii]] += values[ii]
		maps_count[inds[ii]-1][pixels[ii]] += 1
	return v

# Needs cosmology defined from nbodykit, e.g. cosmo=cosmology.Planck15
def z_from_comoving_distance(x,zmax=zmax):
	zgrid = numpy.logspace(-8, 20, 1024)
	zgrid = numpy.concatenate([[0.], zgrid])
	rgrid = cosmo.comoving_distance(zgrid)
  return scipy.interpolate.interp1d(rgrid, zgrid)(x)

# this function takes the density field in fourier space and multiplies with 1/k
def get_x(x, v): return x[0] # get x on the mesh by
def get_y(x, v): return x[1] # get x on the mesh by
def get_z(x, v): return x[2] # get x on the mesh by

# this function takes the density field in fourier space and multiplies with 1/k
def get_vp(k, v):
	kk = sum(ki ** 2 for ki in k) # k^2 on the mesh
	kk[kk == 0] = 1
	return v / np.sqrt(kk) # divide the mesh by k

def get_vx(k, v):
	return 1j * k[0] / k.normp(zeromode=1) ** 0.5 * v

def get_vy(k, v):
	return 1j * k[1] / k.normp(zeromode=1) ** 0.5 * v

def get_vz(k, v):
	return 1j * k[2] / k.normp(zeromode=1) ** 0.5 * v

# this function takes the field (overdensity + 1) in real space and multiplies with aHf
def get_overdensity(x, v):
	return (v-1.)

# this function takes the field in real space and multiplies with aHf
def get_coefficent(x, v):
	x0 = x[0]; y0 = x[1]; z0 = x[2]
	r = np.sqrt(x0**2.+y0**2.+z0**2.)
	tot_l = len(r)*len(r)
	rn = r.reshape(tot_l)
	values = v.reshape(tot_l)
	redshift = zofchi(rn)
	f = cosmo.scale_independent_growth_rate(redshift)
	ah = (100.*cosmo.efunc(redshift))/(1.+redshift)
	faH = (f * ah).reshape(len(r),len(r))
	return v*faH

# v_r = sin(theta)*cos(phi)*v_x + sin(theta)*sin(phi)*v_y + cos(theta)*v_z
def get_v_r_vx(x, v):
	x0 = x[0] + 1e-11
	y0 = x[1]
	z0 = x[2]
	r = np.sqrt(x0**2.+y0**2.+z0**2.)
	# debugging
	x_ = x0*1.+y0*0.+z0*0.
  y_ = x0*0.+y0*0.+z0*0.
	z_ = x0*0.+y0*0.+z0*1.
	p = np.arctan(y_/x_)
	t = np.arccos(z_/r)
	if x0[0,0] > 0:
		p = np.pi - p
	return v * np.sin(t) * np.cos(p)

def get_v_r_vy(x, v):
	x0 = x[0] + 1e-11
	y0 = x[1]
	z0 = x[2]
	r = np.sqrt(x0**2.+y0**2.+z0**2.)
	x_ = x0*1.+y0*0.+z0*0.
  y_ = x0*0.+y0*1.+z0*0.
	z_ = x0*0.+y0*0.+z0*1.
	p = np.arctan(y_/x_)
	t = np.arccos(z_/r)
	if x0[0,0] > 0:
		p = np.pi - p
	return v * np.sin(t) * np.sin(p)

def get_v_r_vz(x, v):
	x0 = x[0] + 1e-11
	y0 = x[1]
	z0 = x[2]
	r = np.sqrt(x0**2.+y0**2.+z0**2.)
	x_ = x0*1.+y0*0.+z0*0.
	y_ = x0*0.+y0*1.+z0*0.
	z_ = x0*0.+y0*0.+z0*1.
	p = np.arctan(y_/x_)
	t = np.arccos(z_/r)
	if x0[0,0] > 0:
		p = np.pi - p
	return - v * np.cos(t)

# v_theta = cos(theta)*cos(phi)*v_x + cos(theta)*sin(phi)*v_y - sin(theta)*v_z
def get_v_theta_vx(x, v):
	x0 = x[0] + 1e-11
	y0 = x[1]
	z0 = x[2]
	r = np.sqrt(x0**2.+y0**2.+z0**2.)
	x_ = x0*1.+y0*0.+z0*0.
	y_ = x0*0.+y0*1.+z0*0.
	z_ = x0*0.+y0*0.+z0*1.
	p = np.arctan(y_/x_)
	t = np.arccos(z_/r)
	if x0[0,0] > 0:
		p = np.pi - p
	return v * np.cos(t) * np.cos(p)

def get_v_theta_vy(x, v):
	x0 = x[0] + 1e-11
	y0 = x[1]
	z0 = x[2]
	r = np.sqrt(x0**2.+y0**2.+z0**2.)
	x_ = x0*1.+y0*0.+z0*0.
  y_ = x0*0.+y0*1.+z0*0.
	z_ = x0*0.+y0*0.+z0*1.
	p = np.arctan(y_/x_)
	t = np.arccos(z_/r)
	if x0[0,0] > 0:
		p = np.pi - p
	return v * np.cos(t) * np.sin(p)

def get_v_theta_vz(x, v):
	x0 = x[0] + 1e-11
	y0 = x[1]
	z0 = x[2]
	r = np.sqrt(x0**2.+y0**2.+z0**2.)
	x_ = x0*1.+y0*0.+z0*0.
  y_ = x0*0.+y0*1.+z0*0.
	z_ = x0*0.+y0*0.+z0*1.
	p = np.arctan(y_/x_)
	t = np.arccos(z_/r)
	if x0[0,0] > 0:
		p = np.pi - p
	return - v * np.sin(t)

# v_phi = -sin(theta)*v_x + cos(phi)*v_y
def get_v_phi_vx(x, v):
	x0 = x[0] + 1e-11
	y0 = x[1]
	z0 = x[2]
	r = np.sqrt(x0**2.+y0**2.+z0**2.)
	x_ = x0*1.+y0*0.+z0*0.
	y_ = x0*0.+y0*1.+z0*0.
	z_ = x0*0.+y0*0.+z0*1.
	p = np.arctan(y_/x_)
	t = np.arccos(z_/r)
	if x0[0,0] > 0:
		p = np.pi - p
	return - v * np.sin(p)

def get_v_phi_vy(x, v):
	x0 = x[0] + 1e-11
	y0 = x[1]
	z0 = x[2]
	r = np.sqrt(x0**2.+y0**2.+z0**2.)
	x_ = x0*1.+y0*0.+z0*0.
  y_ = x0*0.+y0*1.+z0*0.
	z_ = x0*0.+y0*0.+z0*1.
	p = np.arctan(y_/x_)
	t = np.arccos(z_/r)
	if x0[0,0] > 0:
		p = np.pi - p
	return v * np.cos(p)
