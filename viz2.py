import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
plt.rcParams['lines.linewidth'] = 0.5
plt.rcParams['figure.figsize'] = (12, 8)  # Set the figure size to 8x6 inches
plt.rcParams['figure.dpi'] = 300  # Set the DPI to 100


outfolder = 'figs/'
try:
    filepath = sys.argv[1]
except IndexError:
    raise RuntimeError('No file passed')

data = pd.read_csv(filepath)
n = len(data.index)
data = data.iloc[int(n*0.8):]

data.columns = [i.strip() for i in data.columns] # remove whitespaces from colnames
data['d'] = data['b'] - data['a']
data['dtheta'] = data['thetaL'] - data['thetaR']

data.plot(subplots=True)
plt.savefig(outfolder + 'chain.png')
plt.clf()
for c in data.columns:
    data[c].hist()
    plt.title(c)
    plt.savefig(outfolder + c + '.png')
    plt.clf()
