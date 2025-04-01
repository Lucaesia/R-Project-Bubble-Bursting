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
jet = cm = plt.get_cmap('Reds')
N=5
cNorm  = colors.Normalize(vmin=0, vmax=2*N)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

jet = cm = plt.get_cmap('Blues')
N=5
cNorm  = colors.Normalize(vmin=0, vmax=2*N)
scalarMap2 = cmx.ScalarMappable(norm=cNorm, cmap=jet)


# sets time inverval between files
dt = 0.01


#variable for tracking number of files
total_points = 1
scale=2
fig, ax = plt.subplots(5,5,figsize=(5*scale, 5*scale))
# set axes limits manually because Collections do not take part in autoscaling
j=-1
for R in ["0.002","0.0016","0.0012","0.0008","0.0002"]:
  j+=1
  total_points = 1
  i=-1
  for string in ["0.01","0.10","0.20","0.30","0.49"]:
      i += 1
      # finds file name and exits loop once finished
      #R = i
      Level = 10
      IC = "IC"
      t = 10.0
      #string = "%.2f" % float(total_points*0.1)
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
      
      line_collection = mc.LineCollection(lines,edgecolors="red")

      # plots and saves graphs
      ###################################
      # finds file name and exits loop once finished
      #R = "0.002"
      Level = 10
      IC = "NN"
      t = 10.0
      #string = "%.2f" % float(total_points*0.1)
      file_name = "../DriverCode/Water-"+IC+"-R"+R+"-Level"+str(Level)+"/Interfaces/interfaceDrop-"+string+".dat"
      try:
          data = np.loadtxt(file_name, delimiter=' ', usecols=(0,1), unpack=False)
      except FileNotFoundError:
          print("file not found"+file_name)
          break
      # swaps collums and then reagranges into line segments
      data[:, [0, 1]] = data[:, [1, 0]]
      data_mirrored = np.concatenate((data,data*[-1,1]))

      lines = data_mirrored.reshape((len(data_mirrored)//2,2,2))

      colorVal = scalarMap2.to_rgba(total_points)
      # create a LineCollection with the half-circles
      
      line_collection2 = mc.LineCollection(lines,edgecolors="blue")



      ################################
      


      
      ax[j,i].add_collection(line_collection)
      ax[j,i].add_collection(line_collection2)
      ax[j,i].set_xlabel("$x$",fontsize=20)
      ax[j,i].set_ylabel("$z$",fontsize=20)
      ax[j,i].tick_params(axis='both', which='major', labelsize=0,color='w')
      ax[j,i].set_xticks(np.arange(-2, 2, 0.5), minor=True)   # set minor ticks on x-axis
      ax[j,i].set_yticks(np.arange(-2, 2, 0.5), minor=True)   # set minor ticks on y-axis
      ax[j,i].tick_params(which='minor', length=0)       # remove minor tick lines
      ax[j,i].tick_params(which='major', length=0)       # remove minor tick lin
      ax[j,i].set_xlim(-2, 2)
      ax[j,i].set_ylim(-2, 2)
      ax[j,i].grid()                                     # draw grid for major ticks
      ax[j,i].grid(which='minor', alpha=0.3)  
      ax[j,i].set_aspect("equal")  # to make the arcs look circular

      string2 = "%03d" % total_points
      print(string2)
      
      total_points+=2
      if total_points>= N*2:
          break
    

plt.tight_layout()
plt.legend()
plt.savefig('combined_comparison.png')
plt.clf() 
plt.close()

