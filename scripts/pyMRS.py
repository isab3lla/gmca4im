# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 09:40:14 2016

@author: jbobin
"""

## Multi-Resolution on the Sphere: MRS

import numpy as np
import healpy as hpy
from copy import deepcopy

#CXX_LOCATION  = "/Users/jbobin/Documents/Main/exec_healpix"

#
#
#    HEALPIX TOOLBOX IN PYTHON
#
#

def gnside(imag):

    return hpy.get_nside(imag)

#
#  MRS_ALMTRANS : SPHERICAL HARMONICS TRANSFORM (can use the C++ codes or not)
#

def almtrans(map,lmax=None):

    # To be done

    if lmax==None:
        lmax = 3.*hpy.get_nside(map)
        print("lmax = ",lmax)

    alm = hpy.sphtfunc.map2alm(map,lmax=lmax)

    tab = alm2tab(alm,lmax)

    return tab

#
#  MRS_ALMREC : SPHERICAL HARMONICS RECONSTRUCTION (can use the C++ codes or not)
#

def almrec(tab,nside=512):

    alm = tab2alm(tab)

    map = hpy.alm2map(alm,nside,verbose=False)

    return map

#
#  ALM PRODUCT
#

def alm_product(tab,filt):

    length=np.size(filt)
    lmax = np.shape(tab)[0]

    if lmax > length:
        print("Filter length is too small")

    for r in range(lmax):
        tab[r,:,:] = filt[r]*tab[r,:,:]

    return tab

#
#  ALM2TAB : FROM VECTOR TO TAB
#

def alm2tab(alm,lmax):

    taille = np.size(alm)
    tab = np.zeros((lmax+1,lmax+1,2))

    for r in range(0,taille):
        l,m = hpy.sphtfunc.Alm.getlm(lmax,r)
        tab[l,m,0] = np.real(alm[r])
        tab[l,m,1] = np.imag(alm[r])

    return tab
#
#  TAB2ALM : FROM TAB TO ALM
#

def tab2alm(tab):

    lmax = np.int(np.shape(tab)[0])-1
    taille = np.int(lmax*(lmax+3)/2)+1
    alm = np.zeros((taille,),dtype=complex)

    for r in range(0,taille):
        l,m = hpy.sphtfunc.Alm.getlm(lmax,r)
        alm[r] = np.complex(tab[l,m,0],tab[l,m,1])

    return alm

#
#  GETBEAM
#

def getbeam(fwhm=5, lmax=512):

    tor = 0.0174533
    F = fwhm / 60. * tor
    l = np.linspace(0,lmax,lmax+1)
    ell = l*(l+1)
    bl = np.exp(-ell*F*F /16./np.log(2.))

    return bl

#
# CONVOLVE
#

def convolve(map,beam,lmax=512):

    alm = almtrans(map,lmax=lmax)

    tab = alm_product(alm,beam)

    m = almrec(tab,nside=hpy.get_nside(map))

    return m

#
# Wavelet filter
#

def spline2(size,l,lc):

    res =np.linspace(0,size,size+1)

    res = 2.0 * l * res / (lc *size)
    tab = (3.0/2.0)*1.0 /12.0 * (( abs(res-2))**3 - 4.0* (abs(res-1))**3 + 6 *(abs(res))**3 - 4.0 *( abs(res+1))**3+(abs(res+2))**3)

    return tab

def compute_h(size,lc):

    tab1 = spline2(size,2.*lc,1)
    tab2 = spline2(size,lc,1)
    h = tab1/(tab2+0.000001)
    h[np.int(size/(2.*lc)):size]=0.

    return h

def compute_g(size,lc):

    tab1 = spline2(size,2.*lc,1)
    tab2 = spline2(size,lc,1)
    g = (tab2-tab1)/(tab2+0.000001)
    g[np.int(size/(2.*lc)):size]=1

    return g

#
#  WTTRANS : WAVELET TRANSFORM ON THE SPHERE
#

def wttrans(map,nscale=4,lmax=128):

    ech = 1

    taille = np.size(map)

    alm = almtrans(map,lmax=lmax)
    alm_temp = deepcopy(alm)

    LScale = deepcopy(map)
    nside = hpy.get_nside(map)

    WT = np.zeros((taille,nscale))

    for j in range(nscale-1):

        h = compute_h(lmax,ech)
        # g = compute_g(lmax,ech) # Needed if the difference is computed in the spherical harmonics domain

        alm_temp = alm_product(alm,h)

        m = almrec(alm_temp,nside=nside)

        HScale = LScale - m
        LScale = m

        WT[:,j] = HScale

        ech = 2*ech

    WT[:,nscale-1] = LScale

    return WT

####

#
#
#

def wttrans_getfilters(nscale=4,lmax=128):

    ech = 1

    filt = np.zeros((lmax+1,nscale))
    f = np.ones((lmax+1,))

    for j in range(nscale-1):

        h = compute_h(lmax,ech)

        filt[:,j] = f - h*f

        f = h*f

        ech = 2*ech

    filt[:,nscale-1] = f

    return filt

#
#  WTTRANS : WAVELET TRANSFORM ON THE SPHERE // SINGLE SCALE
#

def wttrans_singlescale(map,nscale=4,lmax=128,scale=0):

    F = wttrans_getfilters(nscale=nscale,lmax=lmax)

    f = F[:,scale]

    WT = convolve(map,f,lmax=lmax)

    return WT

#
#  WTREC : WAVELET RECONSTRUCTION ON THE SPHERE
#

def wtrec(WT):

    m = np.sum(WT,axis=1)

    return m

#
#  get_one_face
#

def get_one_face(Imag,NumFace=0,nested=False):

    npix=np.shape(Imag)[0]
    nside =  hpy.npix2nside(npix)
    taille_face = npix/12
    cote=np.int(np.sqrt(taille_face))
    if (nested!=True):
       NewIm=hpy.reorder(Imag, r2n = True)
    else :
      NewIm=Imag
    index = np.zeros((cote,cote))
    index=np.array([hpy.xyf2pix(nside,x,y,0,True) for x in range(nside)
                                                         for y in range(nside)])
    Face =np.resize(NewIm[index+taille_face*NumFace],(cote,cote))

    return Face

#
#  put_one_face
#

def put_one_face(Imag,Face,NumFace=0,nested=False):

   npix=np.size(Imag)
   nside =  hpy.npix2nside(npix)
   taille_face = np.int(npix/12.)
   cote=np.int(np.sqrt(taille_face))
   index = np.zeros((cote,cote))
   index=np.array([hpy.xyf2pix(nside,x,y,0,True) for x in range(nside)
                                                         for y in range(nside)])
   Imag[index+taille_face*NumFace]=np.resize(Face,(cote,cote))

   return Imag

#
# get_ideal_beam
#

def getidealbeam(N, lmin=512, lmax=1024, tozero=True):

    bl = np.zeros((N,))
    bl[0:lmin] = 1.
    Np = lmax-lmin - 1
    #x = np.linspace(0,Np-1,Np) / (Np-1)*3.1416/2
    t = spline2( Np, 1, 1)
    if tozero == True:
        bl[lmin:lmax] = t
        bl[lmax::]=0.
    else:
        bl[lmin:lmax] = bl[lmin:lmax] + t * (1. - bl[lmin:lmax])

    return bl

#
# resize
#

def resize(imag,nside_out=128,alm=False,order_in='NESTED'):

    if alm == False:

        imout = hpy.ud_grade(imag,nside_out,order_in=order_in)

    if alm == True:

        nside = hpy.get_nside(imag)
        lmax=np.min([3*nside,4000])
        alm = almtrans(imag,lmax=lmax)
        lmmax = np.min([3*nside_out,4000])
        lmmin = np.int(0.95*lmmax)
        #print(lmmax)
        if nside_out < nside:
            bl = getidealbeam(lmax+1, lmin=lmmin, lmax=lmmax, tozero=True)
            tab = alm_product(alm,bl)
            imout = almrec(tab,nside=nside_out)
        else:
            imout = almrec(alm,nside=nside_out)


    return imout

#
# tab2nest
#

def tab2nest(nside,xin):

    i = np.mod(xin,(nside/2))
    j = np.fix(xin / (nside/2))
    pas = 1
    out = 0
    nside2 = nside

    while (nside2 > 2):

        out = out +np.mod(i , 2)*pas + 2 *pas * np.mod( j , 2)
        i = np.fix(i/2)
        j= np.fix(j/2)
        nside2 = nside2/2
        pas = pas*4

    return out

#
#  TVS
#

def tvs(imag,  unit='', xsize=800, title='', min=None, max=None, sub=None, nest=False):

    hpy.mollview(imag, fig=None, unit=unit, xsize=xsize, title=title, min=min, max=max, sub=None,nest=nest)

#
#   Test
#

def test():

    LOC = "/Users/jbobin/Documents/Python/LENA_DEVL/Database/Planck2D/Spherical/"

    imag = hpy.read_map(LOC+"wmap_imap_r10_yr9_W4_v5.fits")


def get_all_faces(Imag,nested=False):
   npix=np.shape(Imag)[0]
   nside =  hpy.npix2nside(npix)
   taille_face = np.int(npix/12)
   cote=np.int(np.sqrt(taille_face))
   CubeFace = np.zeros((12,cote,cote))
   if (nested!=True):
      NewIm=hpy.reorder(Imag, r2n = True)
   else :
      NewIm=Imag
   index = np.zeros((cote,cote))
   index=np.array([hpy.xyf2pix(nside,x,y,0,True) for x in range(nside)
                                                         for y in range(nside)])
   for face in range(12):
      #print("Process Face {0}".format(face))
      CubeFace[face,:,:] =np.resize(NewIm[index+taille_face*face],(cote,cote))
      #plt.figure(),imshow(numpy.log10(1+CubeFace[face,:,:]*1e6))
      #plt.title("face {0}".format(face)),plt.colorbar()
   return CubeFace

def put_all_faces( CubeFace,nested=False):
   npix=np.size(CubeFace)
   nside =  hpy.npix2nside(npix)
   taille_face = np.int(npix/12)
   cote=np.int(np.sqrt(taille_face))
   Imag = np.zeros((npix))
   index = np.zeros((cote,cote))
   index=np.array([hpy.xyf2pix(nside,x,y,0,True) for x in range(nside)
                                                         for y in range(nside)])
   for face in range(12):
      #print("Process Face {0}".format(face))
      Imag[index+taille_face*face]=np.resize(CubeFace[face,:,:],(cote*cote))

   if (nested!=True):
      NewIm=hpy.reorder(Imag, n2r = True)
   else:
      NewIm=Imag
   return NewIm

#
# MAP ROTATION
#


def map_rotation(x,theta=22.5,phi=22.5,inv=None,PutToZero=False):

    nside =gnside(x)
    npix = hpy.nside2npix(nside)
    pix = np.arange(npix)

    t,p = hpy.pix2ang(nside,pix) #theta, phi

    r = hpy.Rotator(deg=True,rot=[theta, phi],inv=inv)

    map_rot = -1e10*np.ones(npix)

    trot, prot = r(t,p)
    u = hpy.ang2pix(nside=nside,phi = prot,theta=trot)
    map_rot[u] = x

    ind = np.where(map_rot > -1e9)[0]
    jnd = np.where(map_rot < -1e9)[0]

    if PutToZero:

        map_rot[jnd] = np.mean(map_rot)

    else:

        for i in jnd:
            d = np.sqrt((t[ind] - t[i])**2 + (p[ind] - p[i])**2)
            I = np.argsort(d)
            map_rot[i] = np.mean(map_rot[ind[I[0:10]]])

    return map_rot

def multi_map_rotation(x,theta=22.5,phi=22.5,inv=None,PutToZero=True):

    nside =gnside(x[:,0])
    npix = hpy.nside2npix(nside)
    pix = np.arange(npix)
    n_channels = np.shape(x)[1]

    t,p = hpy.pix2ang(nside,pix) #theta, phi

    r = hpy.Rotator(deg=True,rot=[theta, phi],inv=inv)

    map_rot = -1e10*np.ones((npix,n_channels))

    # for i in pix:
    #     trot, prot = r(t[i],p[i])
    #     u = hpy.ang2pix(nside=nside,phi = prot,theta=trot) # Should be computed offline for different lmax
    #     map_rot[u,:] = x[i,:] #this being the rright way round may need double-checking

    trot, prot = r(t,p)
    u = hpy.ang2pix(nside=nside,phi = prot,theta=trot)
    map_rot[u,:] = x

    ind = np.where(map_rot[:,0] > -1e9)[0]
    jnd = np.where(map_rot[:,0] < -1e9)[0]

    if PutToZero:

        map_rot[jnd,:] = 0.

    else:

        for i in jnd: # Could be somehow made parallel
            d = np.sqrt((t[ind] - t[i])**2 + (p[ind] - p[i])**2)
            I = np.argsort(d)
            map_rot[i,:] = np.mean(map_rot[ind[I[0:10]],:])

    return map_rot

# Get specific regions of the sky

def get_region(lst_id=None,name=None):

    if not hasattr(lst_id, "__len__"):
        lst_id=[lst_id]

    ref_region = ["Polaris","Orion","Pipe","Ophiuchus","Taurus","RCrA",
      "Chamaeleon-South","Pyxis", "Aquila","Auriga","RCrA-Tail","Hercules",
      "Libra","Chamaeleon-Musca","Aquila-Rift","Ara","Pisces", "Microscopium",
      "Triangulum","Perseus","Pavo","Galactic-Center"]
    lc=np.array([120.,211.,0.,354.,173.,10.,315.0,240.0,42.,145.,25.,40.,350.,
         300., 18.,336.,133.,15.,325.,143.,336.,0.])
    bc=np.array([27.,-16.,4.5,15.,-15.,-22.,-22.,12.,-15.,0.,-22.,45.,40.,-13.,
         24.,-14.,-37.,-40.,-14,-25,-28,0.])
    dl = np.array([12.,12.,5.5,12.,12.,15.,12.,25.,10.,50.,15.,15.,30.,12.,25.,
         12.,12.,12.,10.,12.,12.,12.])
    db = np.array([12.,12.,5.5,12.,12.,17.,12.,15.,10.,30.,17.,50.,30.,12.,30.,
         12.,12.,12.,7.,12.,12.,12.])

    lst_match=[]
    if(name is not None):
        lst_match=[i for i,k in enumerate(ref_region) if name in k]
        if len(lst_match) ==0:
            clm=difflib.get_close_matches(name,ref_region)
            lst_match=[i for i,k in enumerate(ref_region) if k in clm]

    if len(lst_match) ==0:
        if(lst_id is not None):
            lst_match=[k for k in lst_id if k <len(ref_region) and k >= 0]
        else:
            print("STRING {0} NOT FOUND, USE GALACTIC CENTER INSTEAD".format(name))
            lst_match=[21]
    return [{'Name':ref_region[k],'lc':lc[k],'bc':bc[k],'dl':dl[k],
                                                'db':db[k]} for k in lst_match]


def extract_region(hmap,lst_id=None,region=None,keep_reso=True,minpix=128,
                              iter=0,ortho=False,fwhm_deg=None,reorder_ns=None):

    try:
        dict_match= get_region(lst_id=lst_id,name=region)
        nside=hp.get_nside(hmap)
        nmaps=hp.maptype(hmap)
        print(nmaps)
        #case single map, write it as a list
        if nmaps==0:
            pmap=[hmap]
            nmaps=1
        else:
            pmap=hmap

        if fwhm_deg is not None:
            if hasattr(fwhm_deg, "__len__"):
                assert(len(fwhm_deg) == nmaps),"Not enough fwhm: {0} vs {1}"\
                                    "map(s).".format(len(fwhm_deg),nmaps)
                fwhm_use= fwhm_deg
            else:
                fwhm_use=np.zeros([nmaps])+fwhm_deg
            lmax=np.zeros((nmaps))
            for kmap in range(nmaps):
                fl=hp.sphtfunc.gauss_beam(fwhm_use[kmap]*DtoR,lmax=2048,
                                                                     pol=False)
                llist=(np.where(fl < 1e-6))[0]
                if len(llist)==0:
                    lmax[kmap]=3*nside-1
                else:
                    lmax[kmap]=llist[0]
                if (reorder_ns is not None) and (hp.isnsideok(reorder_ns)):
                    nside=reorder_ns
                sm_map=[smooth_map_reorder_alm(pmap[kmap], fwhm_deg=fwhm_use[kmap],
                              lmax=lmax[kmap],nside_out=nside,iter=iter)\
                              for kmap in range(nmaps)]
        else:
            sm_map=pmap

        patch=[]
        for kreg in range(len(dict_match)):
            rotation=(dict_match[kreg]["lc"],dict_match[kreg]["bc"],0.)
            if keep_reso is True:
                reso_arcmin=hp.nside2resol(nside,arcmin=True)
                nxpix=np.int(max([np.ceil(dict_match[kreg]["dl"]*60./reso_arcmin),minpix]))
                nypix=np.int(max([np.ceil(dict_match[kreg]["db"]*60./reso_arcmin),minpix]))
                print("STAT=",reso_arcmin,nxpix,nypix,rotation)

            else:
                reso=np.min([dict_match[kreg]["dl"]/minpix,
                                                dict_match[kreg]["db"]/minpix])
                reso_arcmin=reso*60.
                nxpix=np.int(dict_match[kreg]["dl"]/reso)
                nypix=np.int(dict_match[kreg]["db"]/reso)
            a=plt.figure()
            nsub=np.int(np.max((np.ceil(np.sqrt(nmaps)),1)))
            patchreg = [np.array(hp.visufunc.gnomview(map=sm_map[kmap],rot=rotation,fig=a,
               coord='G', xsize=nxpix,ysize=nypix, reso=reso_arcmin,gal_cut=0,
               title=dict_match[kreg]["Name"]+" map "+str(kmap),flip='astro',
               format='%.3g',cbar=True,hold=False,margins=None,notext=False,
               sub=(nsub,nsub,kmap+1),return_projected_map=True))
                                                      for kmap in range(nmaps)]

            patch.append(patchreg)
    except Exception as inst:
        print("FAILURE: ",inst)

    return patch
