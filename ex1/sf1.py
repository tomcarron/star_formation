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
Author(s): Tom Carron
"""
print("************************************************************************")
print("Theoretical and Computational Star formation SoSe 2022 Assignment sheet 1 - Tom Carron")
print("************************************************************************")

"""
Part 1 : Jeans Analysis

Cloud 1: molecular Hydrogen, T=20K and density 1e-18 g/cm^3

Cloud 2: Neutral Hydrogen, T=100K and density 1e-23 g/cm^3

1: What is the Jeans Radius and the corresponding Jeans mass for each cloud?
"""
T_1=20
mu_1=2.38*const.u
rho_1=1e-15 #kg/m^3

T_2=100
mu_2=1.4*const.u
rho_2=1e-20 #kg/m^3

#Functions to calculate the Jeans radius and Jeans mass and the isothermal sound speed. All parameters in SI units.

#Function to calculate the Isothermal sound speed from the temperature and mean molecular weight.
#Physical constants taken from astropy.constants package in SI units. Returns sound speed in m/s.
def sound_speed(T,mu):
    y=np.sqrt(const.k_B*T/mu*const.m_p)
    return y

#Function to calculate the jeans radius and the jeans mass.
#Takes temperature, mean molecular weight and density as inputs in SI units.
#Returns jeans radius in meters and jeans mass in kg, as a two element list [radius,mass].
def jeans(T,mu,rho):
    c_s=sound_speed(T,mu)
    y1=c_s / const.G * rho
    y2=(4/3)*(np.pi)*(rho)*(y1**3)
    return y1.value, y2.value

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
densities=np.logspace(-30,0,25)

plt.figure(0)
plt.loglog(densities,jeans(T_1,mu_1,densities)[1],label="Cloud 1")
plt.loglog(densities,jeans(T_2,mu_2,densities)[1],label="Cloud 2")
plt.xlabel("rho [kg/m^3]")
plt.ylabel("Jeans mass [kg]")
plt.legend()
plt.savefig("plots/rho_jeans_m.png",dpi=400,bbox_inches="tight")
plt.show()
