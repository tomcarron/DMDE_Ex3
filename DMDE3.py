#Dmde Ex3 script
import numpy as np
import matplotlib.pyplot as plt
import astropy
import astropy.units as u
import astropy.constants as const
from astropy.cosmology import Planck18_arXiv_v2 as cosmo
import pandas as pd

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
    term1=(1+(r/(7*rc)))**-1.6
    term2=((2*(r**2) ) / ((rc**2)+(r**2)) )
    term3=(1.6*r / ((7*rc)+r) )
    M = prefactor * r * term1 * (term2 + term3)
    return M

r=np.linspace(0,100,100) #r in Mpc

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

rho500=rho_c*500

'''
Problem 2: Galaxy CLuster mass
For this problem, use the data found in Mgas Mstell.txt on eCampus.
'''
file=open("Mgas_Mstell.txt","r")
'''
data=file.read().split("\n\n")
data=[i.split('\n') for i in data]
df=pd.DataFrame(data,columns=['Name','M_gas(10^13)', 'err_M_gas(10^13)', 'M_tot(10^14)', 'err_M_tot(10^14)', 'M_st(10^12)', 'err_M_st(10^12)'])
print(df)
'''
lines=file.readlines()
name=[]
Mgas=[]
err_Mgas=[]
Mtot=[]
err_Mtot=[]
Mst=[]
err_Mst=[]
for x in lines[1:]:
    name.append(x.split('\t')[0])
    Mgas.append(float(x.split('\t')[1]))
    err_Mgas.append(float(x.split('\t')[2]))
    Mtot.append(float(x.split('\t')[3]))
    err_Mtot.append(float(x.split('\t')[4]))
    Mst.append(float(x.split('\t')[5]))
    err_Mst.append(float(x.split('\t')[6]))
file.close()
print(name)
print(Mst)

plt.figure(1)
plt.scatter(Mgas,Mtot,s=5,color='red')
plt.xlabel('$M_{gas}$ $10^{13}$')
plt.ylabel('$M_{tot}$ $10^{14}$')
plt.grid()
#plt.title('')
#plt.legend()


plt.figure(2)
plt.scatter(Mgas,Mst,s=5,color='green')
plt.xlabel('$M_{gas}$ $10^{13}$')
plt.ylabel('$M_{*}$ $10^{12}$')
plt.grid()
#plt.title('')
#plt.legend()

'''
Plot the gas fraction fgas and the baryon fraction fb vs. Mgas and comment on the
plots. Note that you will need to calculate the fractions from the provided data.
'''

fgas=(1e13*np.array(Mgas))/(1e14*np.array(Mtot))
fb=(1e12*np.array(Mst)+1e13*np.array(Mgas))/(1e14*np.array(Mtot))

plt.figure(3)
plt.scatter(Mgas,fgas,s=5,color='red')
plt.xlabel('$M_{gas}$ ')
plt.ylabel('$f_{gas}$ ')
plt.grid()
#plt.title('')
#plt.legend()


plt.figure(4)
plt.scatter(Mgas,fb,s=5,color='green')
plt.xlabel('$M_{gas}$')
plt.ylabel('$f_{b}$ ')
plt.grid()
#plt.title('')
#plt.legend()


plt.show()
