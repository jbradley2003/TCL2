import os
import pickle
import matplotlib.pyplot as plt
import numpy
import time
from functions import *
from methyl import *

def plot(x,y,xlab,ylab,grid=True):
    plt.plot(x,y)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    if grid:
        plt.grid()
    plt.show()


V = np.linspace(1,500,500)
m = np.array([i for i in range(1,200)])
H = []
B = []
W = []
T = []

def exp_plat(x, a, b):
    return a*(1-np.exp(-b*x))

def fit_exp_plat(x,y):
    popt, pcov = curve_fit(exp_plat, x, y)
    return popt[0],popt[1]

# for s in m:
#     e = eigvals(s,50)
#     B.append(len(e))

for v in V:
    # start = time.perf_counter()
    e = eigvals(40,v)
    t = [i for i in e if i < v]
    u = [i for i in e if i > v]
    if len(u) == 0:
        print('small m')
    n = len(t)
    # end = time.perf_counter()
    # dt = end - start
    # T.append(dt*1e3)
    # B.append(m*2+1)
    W.append(n)

a,b = fit_exp_plat(V,W)
plt.plot(V,exp_plat(V,a,b),color='red',linestyle='--',label=r'$N(V_3) = A(1-e^{-BV_3})$'
         '\n'
         r'$A=$'f'{round(a,3)}, 'r'$B=$'f'{round(b,3)}')
plt.plot(V, W)

# # plt.plot(V,[80]*1000,linestyle='--',color='red')
# plt.xlabel(r'$V_3$ (meV)')
plt.ylabel('Number of Tunneling States')
plt.legend()
plt.title(r'Tunneling Ensemble Size ($m=40$)')
plt.grid()
plt.savefig('conv',dpi=30)
plt.show()
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

