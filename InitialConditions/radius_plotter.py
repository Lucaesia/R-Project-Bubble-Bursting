import numpy as np # type: ignore

from matplotlib import pyplot as plt # type: ignore

array = np.loadtxt("radius list.dat")

scale=1.2
fig, ax = plt.subplots(figsize=(6.4*scale, 3.6*scale))

#ax.plot(np.arange(1.6,31.6,0.2),array*0.3848)
print(np.size(np.arange(16,35.8,0.2)))
print(np.stack((np.arange(16,36,0.2),array*0.3848),axis=1))
#ax.set_xlim([-3,3])
#ax.set_ylim([-0.8,0.7])
#ax.set_aspect('equal')
#np.savetxt("output.txt",[sol2.y[0,:50],sol2.t[:50]])
#ax.plot(men_sol.t, men_sol.y[0,:])
ax.set_ylabel('Major Radius of bubble', fontsize=12)
ax.set_xlabel('θ', rotation=0, fontsize=12)
plt.legend()

plt.savefig('radius.png')
plt.clf()
plt.close()  