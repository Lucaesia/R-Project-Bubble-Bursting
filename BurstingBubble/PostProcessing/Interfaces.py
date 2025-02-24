import numpy as np # type: ignore
from matplotlib import pyplot as plt # type: ignore
from matplotlib import collections  as mc # type: ignore


# sets time inverval between files
dt = 0.1


#variable for tracking number of files
total_points = 1
while True:
    # finds file name and exits loop once finished
    R = 0.0071
    Level = 8
    t = 10.0
    string = "%.1f" % float(total_points*0.1)
    file_name = "../DriverCode/Water-R0.0031-Level8/Interfaces/interfaceDrop-"+string+".dat"
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


    # create a LineCollection with the half-circles
    line_collection = mc.LineCollection(lines)

    # plots and saves graphs
    scale=1
    fig, ax = plt.subplots(figsize=(6.4*scale, 3.2*scale))
    # set axes limits manually because Collections do not take part in autoscaling
    ax.set_xlim(-5, 5)
    ax.set_ylim(-5, 5)
    ax.set_aspect("equal")  # to make the arcs look circular



    ax.add_collection(line_collection)

    string2 = "%03d" % total_points
    print(string2)
    plt.savefig('outputs/interface'+string2+'.png')
    plt.clf() 
    plt.close()
   
    total_points+=1


