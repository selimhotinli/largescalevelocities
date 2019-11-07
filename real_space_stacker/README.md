# The moving lens effect

notes by selim.hotinli14@imperial.ac.uk

## The shape of the signal:

A gravitational potential moving with velocity <img src="http://latex.codecogs.com/svg.latex?\mathbf{v}_\perp" border="0"/> transverse to our line-of-sight direction <img src="http://latex.codecogs.com/svg.latex?\hat{\mathbf{n}}" border="0"/> leads to CMB temperature fluctuations given (at lowest order) by 

<img src="http://latex.codecogs.com/svg.latex?\Delta\Theta (\hat{\mathbf{n}}) = \mathbf{v}_\perp \cdot \boldsymbol{\beta}(\chi \hat{\mathbf{n}})" border="0"/> 

where <img src="http://latex.codecogs.com/svg.latex?\Theta = \Delta\,T/T" border="0"/> is the fractional CMB temperature fluctuation, <img src="http://latex.codecogs.com/svg.latex?\chi" border="0"/> is the conformal distance, and  <img src="http://latex.codecogs.com/svg.latex?\boldsymbol{\beta}" border="0"/> is the deflection angle as seen by the lens (see e.g. arXiv:1812.03167).

On small scales, this moving lens signal is a collection of remote dipolar fluctuations centered on galaxies, aligned with the direction of their transverse velocities. The effect will be detected by the Simons Observatory (SO) to high significance.

Moving lens effect is analogous to (non-linera) Sachs-Wolfe effect, and has the form 

<img src="http://latex.codecogs.com/svg.latex?{ \Theta(\hat{\mathbf{n}})=a/c^2\int\textnormal{d}\chi\,\dot{\Phi}(\chi\hat{\mathbf{n}})}\,," border="0"/> 

where <img src="http://latex.codecogs.com/svg.latex?\chi" border="0"/> is the comoving distance and <img src="http://latex.codecogs.com/svg.latex?\Theta=\Delta\,T/T" border="0"/> is the fractional CMB temperature fluctuation and that 

<img src="http://latex.codecogs.com/svg.latex?{\dot{\Phi}(\chi\hat{\mathbf{n}})}=\mathbf{v}\cdot\nabla\Phi\,." border="0"/>

Assuming spherical symmetry, we can write 

<img src="http://latex.codecogs.com/svg.latex?\boldsymbol{\nabla}\Phi=\hat{\mathbf{r}}\!~\Phi'" border="0"/> 

where <img src="http://latex.codecogs.com/svg.latex?\Phi'=\partial\Phi/\partial\,r" border="0"/>. Using an NFW profile for the halos, i.e.

<img src="http://latex.codecogs.com/svg.latex?\Phi(r)=-4\pi\,G\rho_s\,r_s^2\frac{\ln(1+x)}{x}\,," border="0"/>,

We approximate the functional form of the gravitational potential by using the mass profile for a spherically symmetric halo (NFW) with a single parameter that is the mass of the halo. We fix the virial radius as 

<img src="http://latex.codecogs.com/svg.latex?r_\mathrm{vir}:=\left(\frac{G\,M_\odot\,m}{100H^2}\right)^{1/3}\," border="0"/>

where we will use the parameter <img src="http://latex.codecogs.com/svg.latex?m" border="0"/> to refer to halo mass in Solar mass units, i.e. <img src="http://latex.codecogs.com/svg.latex?M_\odot\simeq1.989\times 10^{30}\mathrm{kg}" border="0"/>. 

Another useful quantity is the clustering parameter, 

<img src="http://latex.codecogs.com/svg.latex?c=A\left(\frac{m}{2\times10^{12}h^{-1}}\right)^\alpha(1+z)^\beta\,," border="0"/>

which relates the scale radius, <img src="http://latex.codecogs.com/svg.latex?r_s" border="0"/>, to the virial radius of the halo, via <img src="http://latex.codecogs.com/svg.latex?c=r_\mathrm{vir}/r_s" border="0"/>. For the model parameters <img src="http://latex.codecogs.com/svg.latex?\{A,\alpha,\beta\}" border="0"/>, we use appropriate values from literature, <img src="http://latex.codecogs.com/svg.latex?\{7.85,-0.081,-0.71\}" border="0"/>. 

We use the NFW profile for the density of the halo,

<img src="http://latex.codecogs.com/svg.latex?\rho(r)=\frac{\rho_s}{x(1+x)^2}\,," border="0"/> 

and 

<img src="http://latex.codecogs.com/svg.latex?\Phi(r)=-4\pi\,G\rho_s\,r_s^2\frac{\ln(1+x)}{x}\,," border="0"/>

where <img src="http://latex.codecogs.com/svg.latex?x=r/r_s" border="0"/>, and we can use the equations above to get 

<img src="http://latex.codecogs.com/svg.latex?\rho_s=M\odot\,m\left[-\frac{r_\mathrm{vir}}{r_s+r\mathrm{vir}}-\ln\left(\frac{r_s+r\mathrm{vir}}{r_s}\right)\right]\,." border="0"/>

The partial derivative of the gravitational potential with respect to r can then be written as 

<img src="http://latex.codecogs.com/svg.latex?\Phi'(r)=4\pi\,G\rho_s\,r_s^2\left[\frac{\ln(1+x)}{x^2}-\frac{1}{x(1+x)}\right]\,," border="0"/>

such that the moving lens signal from a halo takes the form <img src="http://latex.codecogs.com/svg.latex?S=A\tau(x)" border="0"/> where 

<img src="http://latex.codecogs.com/svg.latex?A:=-\frac{v}{c}a\frac{8\pi\,G\rho_sr_s^2}{c^2}\cos\theta" border="0"/>

and the radial dependence is found to be 

<img src="http://latex.codecogs.com/svg.latex?\tau(x):=\frac{1}{2x}\left[\Big|\frac{2\mathrm{arcsec(x)}}{\sqrt{x^2-1}}\Big|+\ln\left(\frac{x^2}{4}\right)\right]\,." border="0"/>           (1)

## Real-space stacking and the detection SNR

Following from our calculations above, we can note that the moving lens effect is a very particular pattern (1) aligned with the transverse velocity of the halos. An efficient method to bring the moving lens signal to surface that stacks the CMB patches benefits from well-approximating the velocities the large scale velocities (so 

