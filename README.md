# large scale velocities (under development)

email: selim.hotinli14@imperial.ac.uk 

Building tools for reconstructing and using large-scale velocity modes

Reconstructed modes are very valuable for cosmological inference as their reconstruction noise scales like  <img src="http://latex.codecogs.com/svg.latex?k^2" border="0"/> and provide better signal-to-noise on largest scales. This makes velocity modes very useful for sample-variance cancellation. 

This rep includes various tools and pipelines for reconstruction large-scale velocities from maps of cosmic microwave background (CMB) and galaxy-surveys. 

tools include: (work on progress)

1) 3D velocity reconstruction from halo catalogs
- uses galaxy catalog to compute 3D overdensity field. 
- uses the 3D overdensity field to compute 3D velocity 3-vector. 
- makes *redshift-binned healpix maps* out of 3D velocities

Requires: 
`healpy`
`nbodykit`

2) stacking algoritm for **moving lens effect**
 - uses halo distributions and CMB maps (with moving lens effect).
 - orients CMB patches to be aligned with the velocity field template. 
 - fits a grid with fixed size in <img src="http://latex.codecogs.com/svg.latex?x=x/r_s" border="0"/> where <img src="http://latex.codecogs.com/svg.latex?r_s" border="0"/> is the scale radius. (See related stacker description.)
 - stacks the grids to extract the signal. 

Moving lens signal is a dipolar modulation in the CMB, alighed with the peculiar velocities of gravitational potentials, e.g. arxiv:1812.03167 

<img src="http://latex.codecogs.com/svg.latex?\Delta\Theta(\hat{\mathbf{n}}) = \mathbf{v}_\perp\cdot\boldsymbol{\beta}(\chi \hat{\mathbf{n}})" border="0"/>

where 

<img src="http://latex.codecogs.com/svg.latex?\boldsymbol{\beta}=\int\textnormal{d}\chi\frac{1}{\chi}\nabla\Phi\,," border="0"/>

similar to the lensing deflection angle. (See more notes at the stacker folder.) 

Requires: 
`healpy`

3) moving lens effect map maker: https://github.com/jbmertens/websky_post
 - uses halo catalog to paint maps of the moving-lens effect.

Requires: 
`healpy`

4) spherical harmonic reconstruction of the transverse velocities:  https://github.com/jbmertens/websky_post (work on progress)

Requires: 
`healpy`

