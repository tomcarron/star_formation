import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

"""
Script to import text files outputted by ex2.cpp and plot them
"""
etas=[2,4,6,8,10]
for i in range(len(etas)):
    densities=[]
    file1=open(str(etas[i])+"_densities.txt","r")
    for line in file1.readlines():
        densities.append(float(line))

    positions=[]
    file2=open(str(etas[i])+"_positions.txt","r")
    for line in file2.readlines():
        positions.append(float(line))

    #print(densities)
    #print(positions)

    fig, ax=plt.subplots(1)
    ax.scatter(positions,densities,label="$\eta$="+str(etas[i]))
    ax.set_ylabel("Density")
    ax.set_xlabel("Position")
    ax.legend(loc="lower center")
    plt.savefig("plots/densities_eta"+str(etas[i])+".png",dpi=400,bbox_inches="tight")
    #plt.show()
