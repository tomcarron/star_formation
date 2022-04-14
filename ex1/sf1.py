import numpy as np
import astropy
import astropy.constants as const
import astropy.units as u
import matplotlib.pyplot as plt
import pandas as pd

"""
Theoretical and Computational Star formation SoSe
2022
Assignment sheet 1
Authors: Tom Carron+Timon Danowski
"""
print("************************************************************************")
print("Theoretical and Computational Star formation SoSe 2022 Assignment sheet 1 - Tom Carron + Timon Danowski")
print("************************************************************************")

"""
Part 1 : Jeans Analysis

Cloud 1: molecular Hydrogen, T=20K and density 1e-18 g/cm^3

Cloud 2: Neutral Hydrogen, T=100K and density 1e-23 g/cm^3

1: What is the Jeans Radius and the corresponding Jeans mass for each cloud?
"""
T_1=20
mu_1=2.38
rho_1=1e-15 #kg/m^3

T_2=100
mu_2=1.4
rho_2=1e-20 #kg/m^3

#Critical density of 1e-13 g/cm^3
rho_crit=1e-10 #g/cm^3

#Functions to calculate the Jeans radius and Jeans mass and the isothermal sound speed. All parameters in SI units.

#Function to calculate the Isothermal sound speed from the temperature and mean molecular weight.
#Physical constants taken from astropy.constants package in SI units. Returns sound speed in m/s.
def sound_speed(T,mu):
    y=np.sqrt(const.k_B*T/(mu*const.m_p))
    return y

#Function to calculate the jeans radius and the jeans mass.
#Takes temperature, mean molecular weight and density as inputs in SI units.
#Returns jeans radius in meters and jeans mass in kg, as a two element list [radius,mass].
def jeans(T,mu,rho):
    c_s=sound_speed(T,mu)
    y1=c_s / np.sqrt(const.G * rho)
    y2=(4/3)*(np.pi)*(rho)*(y1**3)
    return y1.value, y2.value

#barotropic eos function. Returns T in K.
def barotropic_eos(T0,rho,rho_crit,gamma):
    y=T0*(1+(rho/rho_crit)**(gamma-1))
    return y
"""
def jeans_mass(T,mu,rho):
    rad=jeans_radius*T,mu,rho)
    y=(4/3)*(np.pi)*(rho)*(rad**3)
    return y
"""

#Calculating the Jeans mass and radius using the above functions and outputting.
cloud1=jeans(T_1,mu_1,rho_1)
cloud2=jeans(T_2,mu_2,rho_2)

print("************************************************************************")
print("Cloud 1 Jeans radius:",cloud1[0],"[m], Jeans mass",cloud1[1],"[kg]")
print("------------------------------------------------------------------------")
print("Cloud 2 Jeans radius:",cloud2[0],"[m], Jeans mass",cloud2[1],"[kg]")
print("************************************************************************")

"""
2. Plot the Jeans mass a function of density for both clouds on a log-log plot. You
may use Python or any other programing language to produce this plot.

Taking a typical range of densities and evaluating the jeans mass for each density and plotting

"""
densities=np.logspace(-25,-1,48)
density_1=np.logspace(-25,-13,24)
density_2=np.logspace(-13,-1,24)
density=np.logspace(-25,-1,48)

plt.figure(0)
plt.loglog(densities,jeans(T_1,mu_1,densities)[1],label="$H_2$")
plt.loglog(densities,jeans(T_2,mu_2,densities)[1],label="$H$")
plt.xlabel("Density [$kg$ $m^{-3}$]")
plt.ylabel("Jeans mass [$kg$]")
plt.legend()
plt.savefig("plots/rho_jeans_m.png",dpi=400,bbox_inches="tight")


"""
3. Describe your findings.
"""
"""
4. Now consider that the gas behaviour changes at a critical volume density, ρCRIT =
10−13 g/cm 3 . Above this density the adiabatic index is 7/5. Continue your plot to
densities above 10−13 g/cm3 and describe the behaviour of the Jeans mass at these
high densities.
"""
#cloud1_high_density=np.array(jeans(T_1,mu_1,density_1)[1]),np.array(jeans(barotropic_eos(T_1,density_2,rho_crit,7/5),mu_1,density_2)[1])
#cloud2_high_density=np.array(jeans(T_2,mu_2,density_1)[1]),np.array(jeans(barotropic_eos(T_2,density_2,rho_crit,7/5),mu_2,density_2)[1])


cloud1_high_density=[]
cloud2_high_density=[]
i=0
while i < (len(density)):
    if density[i]<1e-13:
        cloud1_high_density.append(jeans(T_1,mu_1,density[i])[1])
        cloud2_high_density.append(jeans(T_2,mu_2,density[i])[1])
        i+=1
    else:
        cloud1_high_density.append(jeans(barotropic_eos(T_1,density[i],rho_crit,7/5),mu_1,density[i])[1])
        cloud2_high_density.append(jeans(barotropic_eos(T_2,density[i],rho_crit,7/5),mu_2,density[i])[1])
        i+=1



plt.figure(1)
plt.loglog(densities,jeans(T_1,mu_1,densities)[1],label="$H_2$")
plt.loglog(densities,jeans(T_2,mu_2,densities)[1],label="$H$")
plt.loglog(density,cloud1_high_density,label="$H_2$-modified")
plt.loglog(density,cloud2_high_density,label="$H$-modified")
plt.xlabel("Density [$kg$ $m^{-3}$]")
plt.ylabel("Jeans mass [$kg$]")
plt.legend()
plt.savefig("plots/rho_jeans_m_high_density.png",dpi=400,bbox_inches="tight")


"""
Q2 after derivation: Calculate the free fall time for both clouds and briefly comment...
"""
#Function to calculate the free fall time of a cloud given its density. Return fft in s. takes rho in kg/m^3.
def free_fall_time(rho):
    y=np.sqrt((3*np.pi) / (32*const.G*rho))
    return y.value

print("------------------------------------------------------------------------")
print("Cloud 1 free-fall time",free_fall_time(rho_1),"[s]")
print("------------------------------------------------------------------------")
print("Cloud 2 free-fall time",free_fall_time(rho_2),"[s]")
print("------------------------------------------------------------------------")
plt.show()
