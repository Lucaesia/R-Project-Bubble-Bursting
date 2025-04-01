import numpy as np # type: ignore
from matplotlib import pyplot as plt # type: ignore
from matplotlib import collections  as mc # type: ignore
from matplotlib import colormaps # type: ignore
import matplotlib.colors as colors # type: ignore
import matplotlib.cm as cmx # type: ignore
import mpl_toolkits.mplot3d.art3d as art3d # type: ignore





# sets time inverval between files
dt = 0.1


#variable for tracking number of files
total_points = 1
scale=0.5
fig, ax = plt.subplots(figsize=(6.4*scale, 3.2*scale))
# set axes limits manually because Collections do not take part in autoscaling
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
ax.set_aspect("equal")  # to make the arcs look circular
while True:
    # finds file name and exits loop once finished
    R = "0.0031"
    Level = 11
    t = 10.0
    string = "%.1f" % float(2.6)
    file_name = "../DriverCode/Water-R"+R+"-Level"+str(Level)+"/Interfaces/interfaceDrop-"+string+".dat"
    try:
        data = np.loadtxt(file_name, delimiter=' ', usecols=(0,1), unpack=False)
    except FileNotFoundError:
        print("file not found"+file_name)
        break
    
      
    ### JET VEL ###
    ax = plt.figure().add_subplot(projection='3d')
    ax.view_init(elev=10)
    ax.set_aspect("equal")  # to make the arcs look circular
    ax.set(xlim=(-3, 3), ylim=(-3, 3), zlim=(-2, 2))

    # swaps collums and then reagranges into line segments
    data[:, [0, 1]] = data[:, [1, 0]]
    data_mirrored = data
    print(data_mirrored)
    R = data_mirrored[:,0]
    Z = data_mirrored[:,1]
    jet = cm = plt.get_cmap('jet')
    N=1
    color_R = np.delete(R, np.arange(0, R.size, 2))
    print(color_R)
    cNorm  = colors.Normalize(vmin=0, vmax=N)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    colors_lines = scalarMap.to_rgba(color_R)


    for i in np.linspace(0,2*np.pi,20):
        x = R*np.cos(i)
        y = R*np.sin(i)
        d3_lines = np.c_[x,y,Z]
        lines = d3_lines.reshape((len(data_mirrored)//2,2,3))

        

        line_collection = art3d.Line3DCollection(lines,colors=colors_lines)
        #ax.add_collection(line_collection)
    R_new = R[R<2]
    Z_new = Z[R<2]
    ring_R = R_new#np.delete(R_new, np.arange(0, R_new.size, 1))
    ring_Z = Z_new#np.delete(Z_new, np.arange(0, R_new.size, 1))
    print(np.where(ring_Z>0.4))
    arr1inds = ring_R.argsort()
    ring_R = ring_R[arr1inds[::-1]]
    ring_Z = ring_Z[arr1inds[::-1]]

    cNorm  = colors.Normalize(vmin=0, vmax=N)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    colors_lines = scalarMap.to_rgba(ring_R)
    
    for i in range(len(ring_Z)):
        theta = np.linspace(0,np.pi,200)
        ax.plot(ring_R[i]*np.cos(theta),ring_R[i]*np.sin(theta),ring_Z[i],color=colors_lines[i],alpha=0.07)
        
          
          

    #d3_lines = np.c_[data_mirrored,np.zeros(len(data_mirrored))]
    lines = d3_lines.reshape((len(data_mirrored)//2,2,3))


    line_collection = art3d.Line3DCollection(lines)
    

    

    

    



    
    

    string2 = "%03d" % total_points
    print(string2)
    
    
    break
    


#plt.legend()
ax.axis('off')
plt.savefig('3D.png', dpi=600)
plt.clf() 
plt.close()

