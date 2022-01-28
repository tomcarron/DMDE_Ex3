#Dmde Ex3 script
import numpy as np
import matplotlib.pyplot as plt
import astropy
import astropy.units as u
import astropy.constants as const
from astropy.cosmology import Planck18_arXiv_v2 as cosmo
import pandas as pd
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


plt.show()
