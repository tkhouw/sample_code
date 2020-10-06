import numpy as np
import pandas as pd

import fileSearcher as fs

runDir = "/mnt/e/runs/"
runName = "squareLattice/R1/N4900/l-3"

dataDir = runDir+runName+"/data/"
searchResults = fs.find("wealthDistribution*0.txt", dataDir)

print(searchResults)

if not searchResults:
    raise OSError("No wealth distributions found in "+dataDir)

bigDF = pd.DataFrame(columns=["w"])
for i in range(len(searchResults)):
    dist = searchResults[i]
    if not (i % 100):
        print(str(int(100*i/len(searchResults)))+"\% completed")
    df = pd.DataFrame()
    wealths = np.loadtxt(dist)
    df["w"] = wealths
    bigDF = bigDF.append(df, ignore_index=True)

print("1\% quantile:")
print(bigDF.quantile(0.01))
print("10\% quantile:")
print(bigDF.quantile(0.1))
print("90\% quantile:")
print(bigDF.quantile(0.9))
print("99\% quantile:")
print(bigDF.quantile(0.99))
print("99.9\% quantile:")
print(bigDF.quantile(0.999))
print("1 - 1/(100*N) quantile:")
print(bigDF.quantile((1 - 1/(100*4900))))
print("1 - 1/(1,000*N) quantile:")
print(bigDF.quantile((1 - 1/(1000*4900))))
print("1 - 1/(10,000*N) quantile:")
print(bigDF.quantile((1 - 1/(10000*4900))))


