import os
import pickle
import matplotlib.pyplot as plt
import numpy
import time
from functions import *
from methyl import *

V = np.linspace(1,1000,800)
m = np.array([i for i in range(1,1000)])
H = []
B = []
W = []
T = []



    
for v in V:
    # start = time.perf_counter()
    e = eigvals(500,v)
    t = [e for e in eigvals if e < v]
    n = len(t)
    # end = time.perf_counter()
    # dt = end - start
    # T.append(dt*1e3)
    # B.append(m*2+1)
    W.append(n)

plt.plot(W, V)
# plt.plot(V,[80]*1000,linestyle='--',color='red')
plt.xlabel('Potential Barrier')
plt.ylabel('Number of Tunneling States')
plt.title(r'Tunneling Ensemble Size ($m= 500$)')
plt.grid()
# plt.savefig('conv',dpi=30)
# plt.show()
# plt.clf()

# plt.plot(V,100*(groundSplittingApprox(V)/(2*np.pi)-W)/W)
# plt.xlabel('Barrier Height (meV)')
# plt.ylabel('Percent Error')
# plt.title('High Barrier Ground Splitting Approximation Error')
# plt.grid()
# plt.show()


'''
import matplotlib.pyplot as plt

fig, ax1 = plt.subplots()

# Create a secondary y-axis that shares the x-axis
ax2 = ax1.twinx()

# Shift the secondary axis to the left
ax2.spines['left'].set_position(('axes', -0.15))
ax2.spines['left'].set_visible(True)
ax2.yaxis.set_label_position('left')
ax2.yaxis.tick_left()

# Hide default right spines
ax2.spines['right'].set_visible(False)

# Add your plots, custom tick colors, etc.

'''

