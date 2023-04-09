import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import astropy.io.fits as fits

nside = 256
npix = nside**2*12
pix = np.arange(npix)
m = np.zeros(npix)

dec_lo = np.pi/2-np.radians(-64.)
dec_hi = np.pi/2-np.radians(18.)

th, phi = hp.pix2ang(nside,pix)
ind = np.where((th>=dec_hi)&(th<=dec_lo))

m[ind] = 1
m[m==0] = -0.1

m_sm = hp.smoothing(m,fwhm=np.radians(3.))
m_sm[m_sm<0] = 0
m_sm/=np.max(m_sm)

print(m_sm)

#m_sm = np.array([m_sm])

print(m_sm)

print(m_sm.shape)

hp.mollview(m_sm)
plt.savefig('ccat_uniform_coverage_nside%i_201021.png'%nside)
plt.show()
plt.clf()

hp.write_map('ccat_uniform_coverage_nside%i_201021.fits'%nside,m_sm,overwrite=True)

