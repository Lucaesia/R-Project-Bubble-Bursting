import numpy as np # type: ignore
from matplotlib import pyplot as plt # type: ignore
from matplotlib import collections  as mc # type: ignore
from matplotlib import colormaps # type: ignore
import matplotlib.colors as colors # type: ignore
import matplotlib.cm as cmx # type: ignore
import matplotlib # type: ignore
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

#plt.style.use('dark_background')


scale=1.5
fig, ax = plt.subplots(figsize=(6.4*scale, 3.2*scale))

ax.plot([2,1.8,1.6,1.4,1.2,1,0.8,0.6,0.4,0.2],[0.225,0.205,0.2,0.175,0.155,0.13,0.107,0.095,0.07,0.03],label = "Using calculated initial condition")
ax.plot([2,1.8,1.6,1.4,1.2,1,0.8,0.6,0.4,0.2],[0.3,0.27,0.25,0.22,0.195,0.16,0.13,0.1,0.065,0.035],label = "Using simplified initial condition")
    
  

plt.legend(fontsize=17)
ax.tick_params(axis='both', which='major', labelsize=17)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.set_xlabel('Bubble Radius', fontsize=17)
ax.set_ylabel('Time Jet Forms', fontsize=17)
plt.savefig('time until.png')
plt.clf() 
plt.close()

