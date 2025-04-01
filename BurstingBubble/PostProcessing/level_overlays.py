import numpy as np # type: ignore
from matplotlib import pyplot as plt # type: ignore
from matplotlib import collections  as mc # type: ignore


# sets time inverval between files
time = 5.0


#variable for tracking number of files
total_points = 1
j=0
scale=4
fig, ax = plt.subplots(6,6,figsize=(6.4*scale, 3.2*scale))
for time in [0.1,1.0,2.0,4.0,6,8]:
    # finds file name and exits loop once finished
    R = "0.0031"
    IC = ""
    Level = 6
    t = 10.0
    string = "%.1f" % float(time)
    file_name = "../DriverCode/Water"+IC+"-R"+R+"-Level"+str(Level)+"/Interfaces/interfaceDrop-"+string+".dat"
    try:
        level_6 = np.loadtxt(file_name, delimiter=' ', usecols=(0,1), unpack=False)
    except FileNotFoundError:
        print("file not found"+file_name)


    Level = 7
    file_name = "../DriverCode/Water"+IC+"-R"+R+"-Level"+str(Level)+"/Interfaces/interfaceDrop-"+string+".dat"
    try:
        level_7 = np.loadtxt(file_name, delimiter=' ', usecols=(0,1), unpack=False)
    except FileNotFoundError:
        print("file not found"+file_name)


    Level = 8
    file_name = "../DriverCode/Water"+IC+"-R"+R+"-Level"+str(Level)+"/Interfaces/interfaceDrop-"+string+".dat"
    try:
        level_8 = np.loadtxt(file_name, delimiter=' ', usecols=(0,1), unpack=False)
    except FileNotFoundError:
        print("file not found"+file_name)


    Level = 9
    file_name = "../DriverCode/Water"+IC+"-R"+R+"-Level"+str(Level)+"/Interfaces/interfaceDrop-"+string+".dat"
    try:
        level_9 = np.loadtxt(file_name, delimiter=' ', usecols=(0,1), unpack=False)
    except FileNotFoundError:
        print("file not found"+file_name)
        

    Level = 10
    file_name = "../DriverCode/Water"+IC+"-R"+R+"-Level"+str(Level)+"/Interfaces/interfaceDrop-"+string+".dat"
    try:
        level_10 = np.loadtxt(file_name, delimiter=' ', usecols=(0,1), unpack=False)
    except FileNotFoundError:
        print("file not found"+file_name)
      

    Level = 11
    file_name = "../DriverCode/Water"+IC+"-R"+R+"-Level"+str(Level)+"/Interfaces/interfaceDrop-"+string+".dat"
    try:
        level_11 = np.loadtxt(file_name, delimiter=' ', usecols=(0,1), unpack=False)
    except FileNotFoundError:
        print("file not found"+file_name)
        
        
          
    ### JET VEL ###


    # swaps collums and then reagranges into line segments
    level_6[:, [0, 1]] = level_6[:, [1, 0]]
    data_mirrored = np.concatenate((level_6,level_6*[-1,1]))
    lines = data_mirrored.reshape((len(data_mirrored)//2,2,2))
    line_collection_6 = mc.LineCollection(lines,colors="red")
    # swaps collums and then reagranges into line segments
    level_7[:, [0, 1]] = level_7[:, [1, 0]]
    data_mirrored = np.concatenate((level_7,level_7*[-1,1]))
    lines = data_mirrored.reshape((len(data_mirrored)//2,2,2))
    line_collection_7 = mc.LineCollection(lines)
    # swaps collums and then reagranges into line segments
    level_8[:, [0, 1]] = level_8[:, [1, 0]]
    data_mirrored = np.concatenate((level_8,level_8*[-1,1]))
    lines = data_mirrored.reshape((len(data_mirrored)//2,2,2))
    line_collection_8 = mc.LineCollection(lines)
    # swaps collums and then reagranges into line segments
    level_9[:, [0, 1]] = level_9[:, [1, 0]]
    data_mirrored = np.concatenate((level_9,level_9*[-1,1]))
    lines = data_mirrored.reshape((len(data_mirrored)//2,2,2))
    line_collection_9 = mc.LineCollection(lines)
    # swaps collums and then reagranges into line segments
    level_10[:, [0, 1]] = level_10[:, [1, 0]]
    data_mirrored = np.concatenate((level_10,level_10*[-1,1]))
    lines = data_mirrored.reshape((len(data_mirrored)//2,2,2))
    line_collection_10 = mc.LineCollection(lines)
    # swaps collums and then reagranges into line segments
    level_11[:, [0, 1]] = level_11[:, [1, 0]]
    data_mirrored = np.concatenate((level_11,level_11*[-1,1]))
    lines = data_mirrored.reshape((len(data_mirrored)//2,2,2))
    line_collection_11 = mc.LineCollection(lines)
    # plots and saves graphs
    
    # set axes limits manually because Collections do not take part in autoscaling
    for i in range(6):
      ax[i,j].set_xlim(-2, 2)
      ax[i,j].set_ylim(-2, 2)
      ax[i,j].set_aspect("equal")  # to make the arcs look circular



    ax[0,j].add_collection(line_collection_6)
    ax[1,j].add_collection(line_collection_7)
    ax[2,j].add_collection(line_collection_8)
    ax[3,j].add_collection(line_collection_9)
    ax[4,j].add_collection(line_collection_10)
    ax[5,j].add_collection(line_collection_11)
    j+=1
    print(j)


plt.savefig('level_overlay.png')
plt.clf() 
plt.close()
total_points+=1
    


