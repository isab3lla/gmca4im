import numpy as np
import sys
import healpy as hp
import time

import pyMRS as pym
import gmca 

## speed of light ##
c = 3.0e8  # m/s

#########################################################
############   to crop the sim (to speed-up)  ###########
#########################################################

def nu_ch_f(nu_ch_in,dnu_out):
	du_in = abs(nu_ch_in[-1]-nu_ch_in[-2])
	a1 = nu_ch_in[0] - du_in/2; a2 = nu_ch_in[-1] + du_in/2
	M = int((a2-a1)/dnu_out)
	if (dnu_out*M)!=(a2-a1):
		print('just dnu multiples!')
		sys.exit()
	nu_ch_out = np.linspace(a1+dnu_out/2,a2-dnu_out/2,M)	

	return nu_ch_out


def merging_maps(nu_ch_in,nu_ch_out,maps_in,dnu_out):
	
	deltanu_in = abs(nu_ch_in[-1]-nu_ch_in[-2])
	maps_out  = [0] * len(nu_ch_out)  

	deltanu_out = abs(nu_ch_out[-1]-nu_ch_out[-2])
	N = int(deltanu_out/deltanu_in)
	if (deltanu_in*N)!=deltanu_out:
		print('just dnu multiples!')
		sys.exit()		

	for i in range(len(nu_ch_out)):
		maps_out[i] = sum(maps_in[N*i:N*i+N]) / N
		
	return maps_out



#########################################################
##################   SURVEY PARAMETERS  #################
#########################################################

## for beam and noise calculation

def survey_info(survey_specs):

	if survey_specs == 'MK': # MeerKAT like
		dish_diam = 13.5  # m
		T_inst    = 20.0  # K
		f_sky     = 0.1   # Survey area (sky fraction)
		t_obs     = 4000. # hrs, observing time
		Ndishes   = 64.   # number of dishes

	if survey_specs == 'custom': # you choose
		dish_diam = 13.5  # m
		T_inst    = 20.0  # K
		f_sky     = 0.1   # Survey area (sky fraction)
		t_obs     = 4000. # hrs, observing time
		Ndishes   = 64.   # number of dishes

	## initialise a dictionary with the instrument specifications
	specs_dict = {'dish_diam': dish_diam, 'T_inst': T_inst, 
					'f_sky': f_sky, 't_obs': t_obs, 'Ndishes' : Ndishes}

	return specs_dict


#########################################################
##################   useful functions   #################
#########################################################

## from vector to matrix and viceversa
def alm2tab(alm,lmax):

	size = np.size(alm)
	tab  = np.zeros((lmax+1,lmax+1,2))

	for r in range(0,size):
		l,m = hp.sphtfunc.Alm.getlm(lmax,r)
		tab[l,m,0] = np.real(alm[r])
		tab[l,m,1] = np.imag(alm[r])

	return tab

def tab2alm(tab):

	lmax = np.int(np.shape(tab)[0])-1
	taille = np.int(lmax*(lmax+3)/2)+1
	alm = np.zeros((taille,),dtype=complex)

	for r in range(0,taille):
		l,m = hp.sphtfunc.Alm.getlm(lmax,r)
		alm[r] = np.complex(tab[l,m,0],tab[l,m,1])

	return alm


## getting the spherical harmonic coefficients
## from a map
def almtrans(map_in,lmax=None):

	if lmax==None:
		lmax = 3.*hp.get_nside(map_in)
		print("lmax = ",lmax)

	alm = hp.sphtfunc.map2alm(map_in,lmax=lmax)

	tab = alm2tab(alm,lmax)

	return tab


## convolution:
## multiplying the spherical harmonic coefficients
def alm_product(tab,beam_l):
	length=np.size(beam_l)
	lmax = np.shape(tab)[0]

	if lmax > length:
		print("Filter length is too small")

	for r in range(lmax):
		tab[r,:,:] = beam_l[r]*tab[r,:,:]

	return tab


## from a_lm back to map
def almrec(tab,nside):

	alm = tab2alm(tab)
	map_out = hp.alm2map(alm,nside,verbose=False)

	return map_out


def plot_cl(fmap,verbose=False):
	LMAX = 3*hp.get_nside(fmap)
	cl = hp.anafast(fmap, lmax=LMAX)
	ell = np.arange(len(cl))
	y = ell * (ell + 1) * cl/2.0/np.pi

	if verbose:
		print("l (l+1) C_l /(2pi) [mK^2] vs l")

	return ell, y

#########################################################
###################   gaussian beam   ###################
#########################################################

## angle in radians of the FWHM
def theta_FWHM(nu,dish_diam): # nu in MHz, dish_diam in m
	return c*1e-6/nu/float(dish_diam) # rad

## solid angle of beam in steradian 
def Omega_beam(nu,dish_diam): # nu in MHz, dish_diam in m 
	return np.pi/(4.*np.log(2))*theta_FWHM(nu,dish_diam)**2

## how many beams to cover my survey area (fraction of sky)
def N_beams(f_sky,nu,dish_diam): # nu in MHz, dish_diam in m 
	return 4*np.pi*f_sky/Omega_beam(nu,dish_diam)

## Fourier transform of the gaussian beam
def getBeam(theta_FWHM,lmax): # theta_FWHM in radians
	sigma_b = theta_FWHM/np.sqrt(8.*np.log(2.))

	l = np.linspace(0,lmax,lmax+1)
	ell = l*(l+1)

	return np.exp(-ell*sigma_b*sigma_b/2)

## convolving the map with the beam
## outputs the new map
def convolve(map_in,beam_l,lmax):

	alm = almtrans(map_in,lmax=lmax)
	tab = alm_product(alm,beam_l)
	m = almrec(tab,nside=hp.get_nside(map_in))

	return m


#########################################################
###################   thermal noise   ###################
#########################################################

def T_sky(nu): # K
	return 60.*(300./nu)**2.55  # K

def T_rcvr(nu,T_inst): # K
	temp_sky = T_sky(nu)
	return 0.1* temp_sky + T_inst

def T_sys(nu,T_inst): # K
	return T_rcvr(nu,T_inst) + T_sky(nu)

## final sigma in mK 
def sigma_N(nu,dnu,T_inst,f_sky,t_obs,Ndishes,dish_diam):
	t_obs = t_obs * 3600 # hrs to s
	dnu = dnu * 1.e6 # MHz to Hz

	temp_sys = T_sys(nu,T_inst)  # in K
	A = np.sqrt(N_beams(f_sky,nu,dish_diam)/dnu/t_obs/Ndishes)
	
	return temp_sys * A *1e3  # mK

def noise_map(sigma,nside=512):
	npixels = hp.nside2npix(nside)
	m = np.random.normal(0.0, sigma, npixels)
	return m

#########################################################
###################    GMCA running   ###################
#########################################################

# J is the number of WT scale
# spherical or 2D patch?
def wavelet_transform(X,J=3):

	print('\nWavelet transforming the data . . .')
	start_w = time.time()

	X_wt = np.zeros((len(X),len(X[0])*J))
	LMAX  = 3*hp.npix2nside(len(X[0]))

	for r in range(X.shape[0]):
		temp = pym.wttrans(X[r],nscale=J+1,lmax=LMAX)
		X_wt[r,:] = temp[:,0:J].reshape(1,-1) # trash the coarse scale

	end_w = time.time()
	tw = end_w - start_w
	print('. . completed in %.2f minutes\n'%(tw/60))
	del tw,end_w,start_w

	return X_wt





def run_GMCA(X_wt,AInit,n_s,mints,nmax,L0,ColFixed,whitening,epsi):

	# First guess mixing matrix (could be set to None or not provided at all)
	if AInit is None:
		AInit = np.random.rand(len(X_wt),n_s)

	print('\nNow running GMCA . . .')

	if whitening:

		R = X_wt@X_wt.T
		L,U = np.linalg.eig(R)
		## whitening the data
		
		Q = np.diag(1./(L+epsi*np.max(L)))@U.T
		iQ = U@np.diag((L+epsi*np.max(L)))

		if ColFixed is None:
			CL = None
		else: CL = Q@ColFixed

		start_w = time.time()
		Results = gmca.GMCA(Q@X_wt,n=n_s,mints=mints,nmax=nmax,L0=L0,Init=0,AInit=AInit,ColFixed=CL)
		end_w = time.time()

		Ae = iQ@Results["mixmat"]  # estimated mixing matrix

	else:
		start_w = time.time()
		Results = gmca.GMCA(X_wt,n=n_s,mints=mints,nmax=nmax,L0=L0,Init=0,AInit=AInit,ColFixed=ColFixed)
		end_w = time.time()

		Ae = Results["mixmat"]

	tw = end_w - start_w
	print('. . completed in %.2f minutes\n'%(tw/60))

	return Ae

