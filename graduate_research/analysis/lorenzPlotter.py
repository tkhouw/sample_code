import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rc

import fileSearcher as fs

#Plots Lorenz curve of resulting wealth distributions

dataDir = "/mnt/e/runs/"
pretty=True
matchSizes=True
lambdaLegend=False
customLegend=True

plot_list = ["fullyConnected/N4900/learn/lis1"]

path_list = [dataDir + pl for pl in plot_list]

full_names = []
for p in path_list:
    full_names += fs.find("wealthDistribution-0.txt", p)
print(full_names)

w_list = [len(np.loadtxt(f))*np.loadtxt(f)/np.sum(np.loadtxt(f)) for f in full_names]

w_list = [np.sort(w) for w in w_list]

w_list = [w/np.sum(w) for w in w_list]
w_list = [np.cumsum(w) for w in w_list]

if matchSizes:
    biggestSize = np.max([len(w) for w in w_list])

ys = np.array(list(range(len(w_list[0]))))/len(w_list[0])
if not matchSizes:
    plt.plot(ys) #plotting line with slope one
else:
    plt.plot(biggestSize*np.arange(len(ys))/len(ys), ys)
for i in range(len(w_list)):
    print("maxes: "+full_names[i]+"\t", np.flip(w_list[i][-4:]))
    if matchSizes:
        plt.plot((biggestSize/len(w_list[i]))*np.arange(len(w_list[i])), w_list[i], alpha=0.8)
    else:
        plt.plot(w_list[i], alpha=0.8)
legend = ["diagonal line"]
if not pretty:
    for p in plot_list:
        legend.append(p)
elif lambdaLegend:
    for f in range(len(full_names)):
        brokenName = full_names[f].split("/")
        exp = None
        for piece in brokenName:
            if (len(piece.split("l-")) == 2) and (piece.split("l-")[0] == ""):
                try:
                    exp = -1*float(piece.split("l-")[1])
                except:
                    pass
        if exp:
            legend.append("1-\u03BB = 10^"+str(round(exp,5)))
        else:
            legend.append("1-\u03BB = -0.01")
elif customLegend:
    legend += [r"BA $\langle k \rangle = 2$", r"BA $\langle k \rangle = 4$", r"BA $\langle k \rangle = 16$", r"BA $\langle k \rangle = 64$", "fully connected"]
plt.legend(legend, prop={"size":14})
plt.xlabel("Rank of wealth", size=14)
plt.ylabel("Cumulative fraction of total wealth", size=14)
plt.tick_params(labelsize=14)
plt.tight_layout()
plt.savefig("plots/result.png")
