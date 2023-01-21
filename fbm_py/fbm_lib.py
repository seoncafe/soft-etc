#!/usr/bin/env python
#--
#-- Written by K.I. Seon (2022.09.19).
#--
#-- fbm2d and fbm3d give 2D and 3D fractal data, respectively, with mean = 0 and standard deviation = 1.
#-- fbm3d_ISM or fbm3d_lognormal_ISM produce a 3D cube data to mimic the ISM fractal density, described by a lognormal distribution.
#-- fbm3d_ISM is the same as fbm3d_lognormal_ISM.
#-- See Seon (2012, ApJL, 761, L17) and Seon & Draine (2016, ApJ, 833, 201).
#--
#-- Usage: The followings show examples how to use this module.
# from fbm_lib import fbm2d, fbm3d, fbm3d_ISM, fbm3d_lognormal_ISM
# a = fbm2d(128,128)
# a.writeto('a.fits.gz')
# plt.imshow(a.data, origin='lower')
# b = fbm3d(128,128,128)
# b.writeto('b.fits.gz')
# plt.imshow(b.data[0,:,:])
# plt.imshow(np.sum(b.data, axis=0))
# c = fbm3d_ISM(128,128,128, mach=2.0)
# c = fbm3d_lognormal_ISM(128,128,128, mach=2.0)
# c.writeto('c.fits.gz')
# plt.imshow(c.data[:,10,:], origin='lower')
# plt.imshow(np.sum(c.data, axis=1), origin='lower')

class fbm2d:
  def __init__(self,nx=64,ny=64,slope=2.4,k_cut_index=None,output_Ak=False,
               seed=None,gaussian_amplitude=False,dtype='float32'):
     import numpy as np
     from numpy.fft import fftfreq, irfft2

     #--- generate random realization
     if seed != None: np.random.seed(seed=seed)

     phi = np.random.random((ny,nx)) * 2.0 * np.pi
     ang = np.zeros((ny,nx//2+1))
     ang[1:,1:nx//2+1] = phi[1:,1:nx//2+1] - phi[:0:-1,:nx//2-1:-1]
     ang[0, 1:nx//2+1] = phi[0, 1:nx//2+1] - phi[0, :nx//2-1:-1]
     ang[1:,0]         = phi[1:,0]         - phi[:0:-1,0]
     ang[0,0]          = 0.0
     del phi

     dx = 1.0
     dy = 1.0
     ky = fftfreq(ny, dy)
     kx = fftfreq(nx, dx)[:nx//2+1]
     Ak = np.zeros((ny,nx//2+1), dtype='complex')

     if k_cut_index != None:
        k_cut = ky[np.int32(k_cut_index)]
     else:
        k_cut = 0.0

     for j in np.arange(ny):
        kr      = np.sqrt(kx[:]**2 + ky[j]**2)
        w       = np.where(kr > k_cut)[0]
        Ak[j,w] = (np.cos(ang[j,w]) + np.sin(ang[j,w])*1j)*kr[w] ** (-slope/2.0)

        if gaussian_amplitude == True:
           gauss = np.random.normal(0.0,1.0,size=(nx//2+1))
           Ak[j,:] = Ak[j,:] * gauss

     Ak[0,0] = 0.0
     Anorm   = 2.0*np.sum(np.abs(Ak[:,1:nx//2+1])**2) + np.sum(np.abs(Ak[:,0])**2) + np.sum(np.abs(Ak[:,nx//2])**2)
     Ak      = Ak/np.sqrt(Anorm)
     img     = irfft2(Ak, norm='forward')

     #--- float32
     if dtype != 'float64': img  = np.float32(img)
     self.data   = img
     self.k_cut  = k_cut
     self.slope  = slope
     if output_Ak == True: self.Ak = Ak

  def writeto(self,fits_file=None,overwrite=True):
     from astropy.io import fits
     if fits_file != None:
        fits_file = fits_file.replace('.fits.gz','').replace('.fits','')+'.fits.gz'
        hdr          = fits.Header()
        hdr['k_cut'] = (self.k_cut, 'wavenumber cut')
        hdr['slope'] = (self.slope, 'power spectrum slope')
        hdu = fits.PrimaryHDU(self.data, header=hdr)
        hdu.writeto(fits_file,overwrite=overwrite)

class fbm3d:
  def __init__(self,nx=64,ny=64,nz=64,slope=2.8,k_cut_index=None,output_Ak=False,
               seed=None,gaussian_amplitude=False,dtype='float32'):
     import numpy as np
     from numpy.fft import fftfreq, irfftn
  
     #--- generate random realization
     if seed != None: np.random.seed(seed=seed)
  
     phi = np.random.random((nz,ny,nx)) * 2.0 * np.pi
     ang = np.zeros((nz,ny,nx//2+1))
     ang[1:,1:,1:nx//2+1]     = phi[1:,1:,1:nx//2+1]    - phi[:0:-1,:0:-1,:nx//2-1:-1]
     #ang(2:nx/2+1,2:ny,2:nz) = phi(2:nx/2+1,2:ny,2:nz) - phi(nx:nx/2+1:-1,ny:2:-1,nz:2:-1)
     ang[1:,1:,0]             = phi[1:,1:,0]            - phi[:0:-1,:0:-1,0]
     #ang(1,2:ny,2:nz)        = phi(1,2:ny,2:nz)        - phi(1,ny:2:-1,nz:2:-1)
     ang[1:,0,1:nx//2+1]      = phi[1:,0,1:nx//2+1]     - phi[:0:-1,0,:nx//2-1:-1]
     #ang(2:nx/2+1,1,2:nz)    = phi(2:nx/2+1,1,2:nz)    - phi(nx:nx/2+1:-1,1,nz:2:-1)
     ang[0,1:,1:nx//2+1]      = phi[0,1:,1:nx//2+1]     - phi[0,:0:-1,:nx//2-1:-1]
     #ang(2:nx/2+1,2:ny,1)    = phi(2:nx/2+1,2:ny,1)    - phi(nx:nx/2+1:-1,ny:2:-1,1)
     ang[1:,0,0]              = phi[1:,0,0]             - phi[:0:-1,0,0]
     #ang(1,1,2:nz)           = phi(1,1,2:nz)           - phi(1,1,nz:2:-1)
     ang[0,1:,0]              = phi[0,1:,0]             - phi[0,:0:-1,0]
     #ang(1,2:ny,1)           = phi(1,2:ny,1)           - phi(1,ny:2:-1,1)
     ang[0,0,1:nx//2+1]       = phi[0,0,1:nx//2+1]      - phi[0,0,:nx//2-1:-1]
     #ang(2:nx/2+1,1,1)       = phi(2:nx/2+1,1,1)       - phi(nx:nx/2+1:-1,1,1)
     ang[0,0,0]               = 0.0
     del phi
  
     dx = 1.0
     dy = 1.0
     dz = 1.0
     kz = fftfreq(nz, dz)
     ky = fftfreq(ny, dy)
     kx = fftfreq(nx, dx)[:nx//2+1]
     Ak = np.zeros((nz,ny,nx//2+1), dtype='complex')

     if k_cut_index != None:
        k_cut = ky[np.int32(k_cut_index)]
     else:
        k_cut = 0.0
  
     #kyy = np.zeros((nz,ny))
     #kzz = np.zeros((nz,ny))
     #for j in np.arange(ny): kyy[:,j] = ky[j]
     #for k in np.arange(nz): kzz[k,:] = kz[k]
     #del ky, kz

     #for i in np.arange(nx//2+1):
     #   kr = np.sqrt(kx[i]**2 + kyy**2 + kzz**2)
     #   if i == 0: kr[0,0] = 1.0
     #   Ak[:,:,i] = (np.cos(ang[:,:,i]) + np.sin(ang[:,:,i])*1j)*kr ** (-slope/2.0)
     #   if gaussian_amplitude == True:
     #      gauss   = np.random.normal(0.0,1.0,size=(nz,ny))
     #      Ak[:,:,i] = Ak[:,:,i] * gauss
     #del kyy, kzz

     kxx = np.zeros((ny,nx//2+1))
     kyy = np.zeros((ny,nx//2+1))
     for i in np.arange(nx//2+1): kxx[:,i] = kx[i]
     for j in np.arange(ny):      kyy[j,:] = ky[j]
     del kx, ky
     for k in np.arange(nz):
        kr = np.sqrt(kxx**2 + kyy**2 + kz[k]**2)
        #if k == 0: kr[0,0] = 1.0
        w        = np.where(kr > k_cut)
        Ak[k][w] = (np.cos(ang[k][w]) + np.sin(ang[k][w])*1j)*kr[w] ** (-slope/2.0)
        if gaussian_amplitude == True:
           gauss = np.random.normal(0.0,1.0,size=(ny,nx//2+1))
           Ak[k] = Ak[k] * gauss
     del kxx, kyy
  
     Ak[0,0,0] = 0.0
     Anorm     = 2.0*np.sum(np.abs(Ak[:,:,1:nx//2+1])**2) + np.sum(np.abs(Ak[:,:,0])**2) + np.sum(np.abs(Ak[:,:,nx//2])**2)
     Ak        = Ak/np.sqrt(Anorm)
  
     # norm = forward gives the standard deviation of 1 (no normalization).
     img    = irfftn(Ak, norm='forward')
     #print('0, mean, stddev = ', np.mean(img), np.std(img))
     #img    = img/np.std(img)
     #print('1, mean, stddev = ', np.mean(img), np.std(img))
  
     #--- float32
     if dtype != 'float64': img  = np.float32(img)
     self.data  = img
     self.k_cut = k_cut
     self.slope = slope
     if output_Ak == True: self.Ak = Ak

  def writeto(self,fits_file=None,overwrite=True):
     from astropy.io import fits
     if fits_file != None:
        fits_file = fits_file.replace('.fits.gz','').replace('.fits','')+'.fits.gz'
        hdr          = fits.Header()
        hdr['k_cut'] = (self.k_cut, 'wavenumber cut')
        hdr['slope'] = (self.slope, 'power spectrum slope')
        hdr['sigma'] = (self.sigma, 'standard deviation')
        hdu = fits.PrimaryHDU(self.data, header=hdr)
        hdu.writeto(fits_file,overwrite=overwrite)
     
class fbm3d_ISM:
  def __init__(self,nx=64,ny=64,nz=64,mach=1.0,bvalue=0.4,normalize=True,k_cut_index=None,
               seed=None,gaussian_amplitude=False,dtype='float32'):

     import numpy as np
     par = np.array([ 2.841e-01, -9.168e-01, -9.334e-01,  1.221e+00, -2.546e-01,
                      8.173e-01,  5.994e-01,  1.326e+00, -1.125e+00,  2.119e-01,
                      5.019e-02, -1.468e-01, -3.838e-01,  2.970e-01, -5.417e-02,
                     -4.428e-03,  1.186e-02,  3.111e-02, -2.375e-02,  4.329e-03]).reshape(4,5)

     # See Seon (2012, ApJL, 761, L17) and Seon & Draine (2016, ApJ, 833, 201)
     # bvalue = 1/3, 0.4, 1.0 for solenoidal, natural mixing, and compressive modes
     if bvalue < 0.4:
        bvalue = 1.0/3.0
     elif bvalue > 0.4:
        bvalue = 1.0

     sigma_lnrho = np.sqrt(np.log(1.0+(bvalue*mach)**2))

     if bvalue <= 0.4:
        slope_ln = 3.81*mach**(-0.16)
     else:
        slope_ln = 3.81*mach**(-0.16) + 0.6

     bcoeff = np.zeros(4)
     for j in np.arange(4):
        bcoeff[j] = np.sum(par[j,:] * sigma_lnrho**np.arange(5))
     slope = np.sum(bcoeff[:] * slope_ln**np.arange(4))

     img = fbm3d(nx=nx,ny=ny,nz=nz,k_cut_index=k_cut_index,slope=slope,seed=seed,
                 gaussian_amplitude=gaussian_amplitude,dtype=dtype)

     self.mach        = mach
     self.bvalue      = bvalue
     self.sigma_lnrho = sigma_lnrho
     self.slope_ln    = slope_ln
     self.slope       = slope
     self.k_cut       = img.k_cut
     self.data        = np.exp(img.data * sigma_lnrho)
     if normalize == True:
        self.data = self.data/np.mean(self.data)

  def writeto(self,fits_file=None,overwrite=True):
     from astropy.io import fits
     if fits_file != None:
        fits_file = fits_file.replace('.fits.gz','').replace('.fits','')+'.fits.gz'
        hdr             = fits.Header()
        hdr['mach']     = (self.mach,        'Mach number')
        hdr['bvalue']   = (self.bvalue,      'b (1/3=solenoidal,0.4=natural,1.0=compressive)')
        hdr['siglnrho'] = (self.sigma_lnrho, 'sigma_lnrho')
        hdr['slope_ln'] = (self.slope_ln,    'slope_ln')
        hdr['k_cut']    = (self.k_cut,       'wavenumber cut')
        hdr['slope']    = (self.slope,       'power spectrum slope')
        hdu = fits.PrimaryHDU(self.data, header=hdr)
        hdu.writeto(fits_file,overwrite=overwrite)

class fbm3d_lognormal:
  def __init__(self,nx=64,ny=64,nz=64,slope=2.4,sigma=1.0,normalize=True,k_cut_index=None,
               seed=None,gaussian_amplitude=False,dtype='float32'):

     import numpy as np

     img = fbm3d(nx=nx,ny=ny,nz=nz,k_cut_index=k_cut_index,slope=slope,seed=seed,
                 gaussian_amplitude=gaussian_amplitude,dtype=dtype)

     self.slope = slope
     self.sigma = sigma
     self.k_cut = img.k_cut
     self.data  = np.exp(img.data * sigma)
     if normalize == True:
        self.data = self.data/np.mean(self.data)

  def writeto(self,fits_file=None,overwrite=True):
     from astropy.io import fits
     if fits_file != None:
        fits_file = fits_file.replace('.fits.gz','').replace('.fits','')+'.fits.gz'
        hdr          = fits.Header()
        hdr['slope'] = (self.slope, 'power spectrum slope of log(rho)')
        hdr['sigma'] = (self.sigma, 'sigma_log10')
        hdr['k_cut'] = (self.k_cut, 'wavenumber cut')
        hdu = fits.PrimaryHDU(self.data, header=hdr)
        hdu.writeto(fits_file,overwrite=overwrite)

fbm3d_lognormal_ISM = fbm3d_ISM
