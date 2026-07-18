import os
import numpy as np
from methyl import*

K = 0.0862

def Z(V,T):
    m = 50
    beta = 1/(K*T)
    states = [s for s in eigvals(m, V) if s < V]
    z = 0
    for e in states:
        z += np.exp(-beta*e)
    return z, states

def P(V,T):
    beta = 1/(K*T)
    z, states = Z(V,T)
    print(states)
    P_map = {}
    for n, e in enumerate(states):
        boltz_factor = np.exp(-beta*e)/z
        if boltz_factor in P_map.values():
            key = [i for i in P_map.keys() if P_map[i] == boltz_factor][0]
            P_map[key] = 2*boltz_factor/z
        else:
            P_map[n] = np.exp(-beta*e)/z
    return P_map

if __name__ == "__main__":
    V = 100
    T = np.linspace(5,500,10)

    for t in T:
        p = P(V,t)
        print(p.items())
        plt.plot(p.keys(),p.values(),label=f'T={int(t)} K')

    plt.xticks(list(p.keys()))
    plt.legend()
    plt.show()