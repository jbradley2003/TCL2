import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from parseXYZ import *
import os

INERTIA=5.3e-47
HBAR=1.05457e-34
B_J = HBAR**2/(2*INERTIA)

def createH(m, V_k, k=3, B=B_J*6.242e+21):
  states = [i for i in range(-m, m+1)]
  T = np.diag([B*(i**2) for i in states])
  V = np.diag([V_k/2]*(2*m+1))
  for i in range(2*m+1):
    for j in range(2*m+1):
      if states[i] == states[j] + k or states[i] == states[j] - k:
        V[i,j] = -V_k/4
  H = T + V
  return H, states

def eigvals(m, V):
  return np.linalg.eigvalsh(createH(m, V)[0])

def splitting(m, V):
  evals = eigvals(m, V)
  grouped = [evals[i:i+3] for i in range(0, len(evals), 3)]
  split = [abs(g[0] - g[2]) for g in grouped if len(g) == 3]
  return split, grouped

def nthSplitting(n, V, tol=1e-12):
  m_init = int(np.ceil((3*n-1)/2))
  m = m_init+1
  v_init = splitting(m_init, V)[0][n-1]
  v_curr = splitting(m, V)[0][n-1]
  dv = abs(v_curr - v_init)
  while dv > tol:
       m += 1
       v_next = splitting(m, V)[0][n-1]
       dv = abs(v_next - v_curr)
       v_curr = v_next
  return v_curr, m

def groundSplittingApprox(V, I=INERTIA, hbar=HBAR):
  eV_to_J = 1.602176634e-19
  V_3 = V * 1e-3 * eV_to_J
  const = 6*np.sqrt(3/np.pi)
  return const*(((2*(V_3**3))/(I*hbar**2))**(1/4))*np.exp(-(4*np.sqrt(2*I*V_3))/(3*hbar))

def getBarrier(filename):
  kcal_to_meV = 43.364
  csv = np.loadtxt(filename, delimiter=',', skiprows=1)
  x = np.deg2rad(csv[:,1])
  y = csv[:,4]*kcal_to_meV
  return fitPotential(csv),[x,y]

def getBarriers(root):
  bmap = {}
  # This depends on the following directory name convention: 'Me{central_carbon}' --> make this dynamic
  for d in os.scandir(root):
      if os.path.isdir(d.path) and 'Me' in d.name:
        key = int(d.name[2:])
        csv = [f for f in os.listdir(d.path) if f.endswith('.csv')][0]
        filename = f'{root}/{d.name}/{csv}'
        print(filename)
        barrier, points = getBarrier(filename)
        bmap[key] = (barrier, points)
  return bmap

def fitPotential(data):
  kcal_to_meV = 43.364
  x = np.deg2rad(data[:,1])
  y = data[:,4]*kcal_to_meV
  # guess_barrier = (np.max(y) - np.min(y)) / 2
  popt, pcov = curve_fit(potential, x, y)
  return popt[0]

def potential(phi, barrier):
  return barrier*(1 - np.cos(3*phi))/2

def plotPotentials(root, steps=500):
  bmap = getBarriers(root)
  phi_range = np.linspace(0,2*np.pi,steps)
  for k, v in bmap.items():
      barrier, data = v
      plt.scatter(data[0], data[1], color='red')
      plt.plot(phi_range, potential(phi_range, barrier), color='black', 
                label= r'$V_3$ = ' + str(round(barrier,2)))
      plt.xlabel(r'$\phi$ (rad)')
      plt.ylabel(r'$V(\phi)$ (meV)')
      plt.title(f'Rotational Potential for Me{k}')
      plt.legend()
      plt.show()
      plt.clf()

def methylMap(xyz_path, me_path):
  meV_to_kHz = 241798930
  mol, atoms, xyz = xyzToMol(xyz_path)
  methyls, methyl_bonds = findMethyls(mol)
  bmap = getBarriers(me_path)
  vmap = {}
  for central_carbon, values in methyls.items():
    pairs = values[1]
    for pair in pairs:
      if central_carbon in bmap.keys():
        v = bmap[central_carbon][0]
        v_t = nthSplitting(1,v)[0]
      else:
        print(f'ATOM {central_carbon} is not present in bmap!')
      vmap[pair] = -2*meV_to_kHz*v_t/3
  return vmap