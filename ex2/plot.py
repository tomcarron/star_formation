import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

"""
Script to import text files outputted by ex2.cpp and plot them
"""
densities=[]
file1=open("densities.txt","r")
for line in file1.readlines():
    densities.append(float(line))

positions=[]
file2=open("positions.txt","r")
for line in file2.readlines():
    positions.append(float(line))

#print(densities)
#print(positions)

fig, ax=plt.subplots(1)
ax.scatter(positions,densities)
ax.set_ylabel("Density")
ax.set_xlabel("Position")
plt.savefig("plots/densities.png",dpi=400,bbox_inches="tight")
plt.show()
