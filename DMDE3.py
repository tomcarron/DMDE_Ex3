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
Tgas0=3000#*u.eV#.to('K') #keV
print(Tgas0)

def xray_mass_profile(r,rc,Tgas0):
    rc=rc
    Tgas0=Tgas0
    G=const.G.value
    kB=const.k_B.value #.to(u.cgs).value
    mu=1 #Mean molcular weight of ICM
    mp=const.m_p.value
    prefactor=((3000)/(G*mu*mp))
    term1=(1+(r/(7*rc)))**-1.6
    term2=((2*(r**2) ) / ((rc**2)+(r**2)) )
    term3=(1.6*r / ((7*rc)+r) )
    M = prefactor * r * term1 * (term2 + term3)
    return M

r=np.linspace(0,1,1000) #r in Mpc

#print(const.k_B,const.m_p,const.G)

plt.figure(0)
plt.plot(r,xray_mass_profile(r,rc,Tgas0))
plt.xlabel('r [Mpc]')
plt.ylabel('$M_{tot}(<r)$')
plt.legend()
plt.savefig('plots/fig0.png',dpi=400,bbox_inches="tight")
plt.title('')

'''
Determine r500 and M500 numerically. Remember: r500 is the radius within which the
mean mass density of the cluster is 500 times greater than the critical density of the
Universe, and M500 is the mass within this radius.
'''

rho_c=(const.M_sun.to('kg').value*3*(cosmo.H0.value**2)) / (8*np.pi*(const.G.value)) #units converted to give M_sun/Mpc
#rho_c=(3*(cosmo.H0.value**2)) / (8*np.pi*(const.G.value))
rho_500=rho_c*500
print(rho_500)

#r500 is the radius within which the mean density is 500 * rho_c
r2=np.linspace(0,1e-3,1000)
M_arr=xray_mass_profile(r2,rc,Tgas0)
density_i=0
i=1

while density_i<rho_500:
    M_i=M_arr[i]
    r_i=r2[i]
    volume_i=(4/3)*np.pi*r_i**3
    density_i=M_i/volume_i
    #print(density_i)
    i+=1

print('r500=',r2[i],'M500=',M_arr[i],'i=',i)

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
plt.savefig('plots/fig1.png',dpi=400,bbox_inches="tight")
plt.grid()
#plt.title('')
#plt.legend()


plt.figure(2)
plt.scatter(Mgas,Mst,s=5,color='green')
plt.xlabel('$M_{gas}$ $10^{13}$')
plt.ylabel('$M_{*}$ $10^{12}$')
plt.savefig('plots/fig2.png',dpi=400,bbox_inches="tight")
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
plt.xlim(0,16)
plt.grid()
plt.savefig('plots/fig3.png',dpi=400,bbox_inches="tight")
#plt.title('')
#plt.legend()

dummy=np.linspace(-1,16,30)
universal_fb=np.zeros_like(dummy) + 2/15

plt.figure(4)
plt.scatter(Mgas,fb,s=5,color='green')
plt.plot(dummy,universal_fb,'--',color='red')
plt.xlim(0,16)
plt.xlabel('$M_{gas}$')
plt.ylabel('$f_{b}$ ')
plt.savefig('plots/fig4.png',dpi=400,bbox_inches="tight")
plt.grid()
#plt.title('')
#plt.legend()

print(np.sum(fb)/len(fb)) #mean fb
plt.show()
