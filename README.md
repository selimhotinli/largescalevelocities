# largescalevelocities

tools for reconstructing and using large-scale velocity modes

this rep includes various tools and pipelines for reconstruction large-scale velocities from maps of cosmic microwave background (CMB) and galaxy-surveys. 

tools include: (work on progress)

1) velocity map reconstruction from halo catalogs
- uses galaxy catalog to compute 3D overdensity field. 
- uses the 3D overdensity field to compute 3D velocity 3-vector. 
- makes z-binned healpix maps out of 3D velocities

Requires: 
`healpy`
`nbodykit`

2) stacking algoritm for moving lens effect
 - uses halo distributions and CMB maps (with moving lens effect).
 - orients CMB patches to be aligned with the velocity field template. 
 - fits a grid with fixed size in x=x/r_s where r_s is the scale radius.
 - stacks the grids to extract the signal. 

Requires: 
`healpy`

3) moving lens effect map maker: https://github.com/jbmertens/websky_post
 - uses halo catalog to paint maps of the moving-lens effect.

Requires: 
`healpy`

4) spherical harmonic reconstruction of the transverse velocities:  https://github.com/jbmertens/websky_post (work on progress)

Requires: 
`healpy`

