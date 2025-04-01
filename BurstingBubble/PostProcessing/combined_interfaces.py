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

jet = cm = plt.get_cmap('jet')
N=5
cNorm  = colors.Normalize(vmin=0, vmax=2*N)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)


# sets time inverval between files
dt = 0.1


#variable for tracking number of files
total_points = 1
scale=2
fig, ax = plt.subplots(figsize=(6.4*scale, 3.2*scale))
# set axes limits manually because Collections do not take part in autoscaling
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
ax.set_aspect("equal")  # to make the arcs look circular
while True:
    # finds file name and exits loop once finished
    R = "0.002"
    Level = 9
    IC = "IC"
    t = 10.0
    string = "%.1f" % float(total_points*0.1)
    file_name = "../DriverCode/Water-"+IC+"-R"+R+"-Level"+str(Level)+"/Interfaces/interfaceDrop-"+string+".dat"
    try:
        data = np.loadtxt(file_name, delimiter=' ', usecols=(0,1), unpack=False)
    except FileNotFoundError:
        print("file not found"+file_name)
        break
    
      
    ### JET VEL ###


    # swaps collums and then reagranges into line segments
    data[:, [0, 1]] = data[:, [1, 0]]
    data_mirrored = np.concatenate((data,data*[-1,1]))

    lines = data_mirrored.reshape((len(data_mirrored)//2,2,2))

    colorVal = scalarMap.to_rgba(total_points)
    # create a LineCollection with the half-circles
    if total_points == N or total_points== N+1 or total_points == 1 or total_points>= N*2-2:
        line_collection = mc.LineCollection(lines,edgecolors=colorVal,label="$t$ = "+string)
    else:
      line_collection = mc.LineCollection(lines,edgecolors=colorVal)

    # plots and saves graphs
    


    
    ax.add_collection(line_collection)

    string2 = "%03d" % total_points
    print(string2)
    
    total_points+=1
    if total_points>= N*2:
        break
    

ax.set_xlabel("$x$")
ax.set_ylabel("$z$")
plt.legend()
plt.savefig('combined.png')
plt.clf() 
plt.close()

