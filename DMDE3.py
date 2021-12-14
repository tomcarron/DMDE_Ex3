#Dmde Ex3 script
import numpy as np
import matplotlib.pyplot as plt
import astropy
import astropy.units as u
import astropy.constants as const
from astropy.cosmology import Planck18_arXiv_v2 as cosmo

#Using rc = 0.15 Mpc and Tgas (0) = 3 keV, plot Mtot(< r). Describe and explain the
#shape of the plot.

rc=0.15 #Mpc
Tgas0=3 #keV

def xray_mass_profile(r,rc,Tgas0):
    rc=rc
    Tgas0=Tgas0
    G=const.G.value
    kB=const.k_B.value
    mu=1 #Mean molcular weight of ICM
    mp=const.m_p.value
    prefactor=((kB*Tgas0)/(G*mu*mp))
    term1=((2*(r**2) ) / ((rc**2)+(r**2)) )
    term2=(1.6*r / ((7*rc)+r) )
    M = prefactor * r * (term1 + term2)
    return M

r=np.linspace(0,1,100) #r in Mpc

#print(const.k_B,const.m_p,const.G)

plt.figure(0)
plt.plot(r,xray_mass_profile(r,rc,Tgas0))
plt.xlabel('r [Mpc]')
plt.ylabel('$M_{tot}(<r)$')
plt.legend()
plt.title('')

'''
Determine r500 and M500 numerically. Remember: r500 is the radius within which the
mean mass density of the cluster is 500 times greater than the critical density of the
Universe, and M500 is the mass within this radius.
'''

rho_c=(3*(cosmo.H0.value**2)) / (8*np.pi*(const.G))

#r500 is the radius within which the mean density is 500 * rho_c




plt.show()
