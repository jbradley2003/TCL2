#
#
### imports
#
#

import numpy as np
import pycce as pc
import sympy as sp
import pickle
import math
import re
import copy
from itertools import combinations
from scipy.optimize import curve_fit
from methyl import *
import constants as const

#
#
### miscellaneous functions
#
#

def stretched_exp(xdata,T2,q):
    return 1/2*np.exp(-(xdata/T2)**q)

def exp(xdata,T2):
    return 1/2*np.exp(-xdata/T2)

def fit_stretched_exp(data,time_points,p0=None,verbose=True):
    popt, pcov = curve_fit(stretched_exp,time_points,data,p0=p0)
    if verbose:
        print(f'optimized parameters: {popt}')
        print(f'std dev of parameters: {np.sqrt(np.diag(pcov))}')
        print(f'T2 fit: {popt[0]}')
    elif not verbose:
        print(f'T2 fit: {popt[0]}')
    return stretched_exp(time_points,popt[0],popt[1]), popt

def fit_exp(data,time_points,p0=None,verbose=True):
    popt, pcov = curve_fit(exp,time_points,data,p0=p0)
    if verbose:
        print(f'optimized parameters: {popt}')
        print(f'std dev of parameters: {np.sqrt(np.diag(pcov))}')
        print(f'T2 fit: {popt[0]}')
    elif not verbose:
        print(f'T2 fit: {popt[0]}')
    return exp(time_points,popt[0]), popt

def imapPointDipole(spins,highestN):

    mu0 = np.pi*4e-7 # N A-2
    h = 6.626e-34 # J s

    intes = list(spins.keys())
    
    interactions = combinations(intes,2)

    imap = {}
    distMap = {}
    prefs = {}

    for inte_ in interactions:

        i = inte_[0]
        j = inte_[1]
       
        # record the distances that are computed in the dist function and return as needed        

        dist = np.sqrt((spins[i][0][0]-spins[j][0][0])**2 + (spins[i][0][1]-spins[j][0][1])**2 + (spins[i][0][2]-spins[j][0][2])**2)
        r12 = spins[i][0] - spins[j][0]
        
        prefactor = ((mu0*spins[i][1]*spins[j][1]*h**2) / (16*np.pi**3*dist**3))*(1000/(h*(1e-10)**3*(1e-4)**2))

        aTens = np.array([[(1 - (3*r12[0]**2/dist**2)),-3*r12[0]*r12[1]/dist**2,-3*r12[0]*r12[2]/dist**2],
                    [-3*r12[1]*r12[0]/dist**2,(1 - (3*r12[1]**2/dist**2)),-3*r12[1]*r12[2]/dist**2],
                    [-3*r12[2]*r12[0]/dist**2,-3*r12[2]*r12[1]/dist**2,(1 - (3*r12[2]**2/dist**2))]])

        imap[inte_] = prefactor*aTens
        distMap[inte_] = dist
    
    return imap,distMap

def initializeStatistics(filepath,soi,boxLength,bathType,concentration):

    xyz,hf = readOrcaV2(filepath,soi)

    atoms = pc.random_bath('1H', boxLength,
                           density=concentration,
                           density_units='cm-3')
    
    hCount = 0
    fCount = 0
    for key in list(xyz.keys()):
        if key.split()[-1][-1] == 'F':
            fCount += 1
        else:
            hCount += 1

    numMmPairs = math.comb(fCount,2) + math.comb(hCount,2)
    numMsPairs = hCount*len(atoms)
    numSsPairs = math.comb(len(atoms),2)

    return numMmPairs,numMsPairs,numSsPairs

def initializeStatisticsV2(filepath,soi,boxLength,bathType,concentration):

    xyz,hf = readOrcaV2(filepath,soi)
    
    atoms = pc.random_bath('1H', boxLength,
                           density=concentration,
                           density_units='cm-3')
    
    numBathSpins = len(atoms)
    
    numBathPairs = math.comb(numBathSpins,2)
    
    hCount = 0
    fCount = 0
    keyList = []
    for key in list(xyz.keys()):
        if key.split()[-1][-1] == 'F':
            fCount += 1
        else:
            keyList.append(key)
            hCount += 1

    return numBathSpins,numBathPairs,keyList

def addRandomBath(boxLength,bathType,concentration,spins,imap,molecularAtoms):

    # boxLength is the edge length of a cubic box in angstroms - good value is 30
    # bathType is the bath spin type - 'H' for hydrogen or 'D' for deuterium, nothing else will work right now
    # concentration is the spin concentration in cm^-3 - good values are around 2e22 for a dense bath, could do 1e20 for deuterated proton bath
    
    cpimap = copy.deepcopy(imap)
    cpspins = copy.deepcopy(spins)

    if bathType == 'H':
    
        atoms = pc.random_bath('1H', boxLength, 
                           density=concentration, 
                           density_units='cm-3')
        
    else:
        raise ValueError('only implemented for hydrogen bath currently')
    
    highestN = re.findall(r'\d+', list(spins.keys())[-1])[0]
    
    extraHN = int(highestN) + 1
    
    avgdist = 0
    
    for arr in atoms.xyz:
        avgdist += vec_dist(spins['e'][0],arr)            
        
    avgdist /= len(atoms.xyz)
    
    newavgdist = 0
    
    for xyzarr in atoms.xyz:
        
        xyzarr = spins['e'][0] - xyzarr # correcting for different distances from the spin density on the molecules
        
        newavgdist += vec_dist(spins['e'][0],xyzarr)       

        # this part rejects a bath spin placement if it is too close to any of the molecular atoms

        truthArr = []
        for molecularAtom in molecularAtoms:
            if vec_dist(xyzarr,molecularAtom) <= 1.0:
                truthArr.append(False)
            elif vec_dist(xyzarr,molecularAtom) > 1.0:
                truthArr.append(True)

        if False in truthArr:
            continue
        elif False not in truthArr:
        
            cpspins[f'{extraHN}{bathType}'] = [xyzarr,const.gyro_dict[f'{bathType}']]
            extraHN += 1
    
    newavgdist /= len(atoms.xyz)
    
    imapBath,distMap = imapPointDipole(cpspins,highestN)
    
    for key in list(imapBath.keys()):

        i = key[0]
        j = key[1]
    
        if i == 'e':
            
            jCheck = int(re.findall(r'\d+', j)[0])
            
            if jCheck > int(highestN): 
                
                # this avoids overwriting the hyperfine couplings generated from ORCA
                # and adds all new hyperfine couplings for the random bath
            
                cpimap[(f'{i}',f'{j}')] = imapBath[(i,j)]
                
        elif i != 'e': 
            
            # this adds all new nuclear spin - nuclear spin interactions
            # and keeps old ones, as they're all generated from the point dipole approximation anyways
            
            cpimap[(f'{i}',f'{j}')] = imapBath[(i,j)]
    
    return cpimap,int(highestN),distMap

def vec_dist(v1,v2):
    return np.sqrt((v1[0]-v2[0])**2+(v1[1]-v2[1])**2+(v1[2]-v2[2])**2)

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def vec_magnitude(vec):
    return np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)

def TraceDistance(p1,p2):

    T = 1/2*(np.trace(sqrtm((np.conjugate(p1-p2).T)@(p1-p2))))

    return T

def Fidelity(p1,p2):

    F = (np.trace(sqrtm(sqrtm(p1)@p2@sqrtm(p1))))**2

    return F

def seperate_string_number(string):
    previous_character = string[0]
    groups = []
    newword = string[0]
    for x, i in enumerate(string[1:]):
        if i.isalpha() and previous_character.isalpha():
            newword += i
        elif i.isnumeric() and previous_character.isnumeric():
            newword += i
        else:
            groups.append(newword)
            newword = i

        previous_character = i

        if x == len(string) - 2:
            groups.append(newword)
            newword = ''
    return groups

#
#
### TCL Equation Functions
#
#

# Added J12
A1,A2,b12,J12,wn,wh = sp.symbols('A_1 A_2 b_{12} J_{12} \omega_n \omega_h',real=True)
t = sp.symbols('t',real=True,positive=True)

# Need to make correlation functions in sympy

with open(f'{os.getcwd()}/equations/I1o2I1o2WtTCL2.txt','rb') as f:
    I1o2I1o2WtTCL2 = pickle.load(f)

with open(f'{os.getcwd()}/equations/I1o2I1o2WtTCL2J.txt','rb') as f:
    I1o2I1o2WtTCL2J = pickle.load(f)

with open(f'{os.getcwd()}/equations/I1o2I1o2WtTCL4.txt','rb') as f:
    I1o2I1o2WtTCL4 = pickle.load(f)

# Methyl J-coupling
modulationDepthTCL2J = 4*((b12-2*J12)**2*(A1-A2)**2)/((b12-2*J12)**2 + (A1-A2)**2)**2
frequencyTCL2J = sp.Rational(1,4)*sp.sqrt(A1**2-2*A1*A2+A2**2+(b12-2*J12)**2)

modulationDepthTCL2LambdaJ = sp.lambdify([A1,A2,b12,J12],modulationDepthTCL2J,'numpy')
frequencyTCL2LambdaJ = sp.lambdify([A1,A2,b12,J12],frequencyTCL2J,'numpy')

I1o2I1o2WtTCL2LambdaJ = sp.lambdify([A1,A2,b12,J12,t],I1o2I1o2WtTCL2J,'numpy')

# Original TCL2
modulationDepthTCL2 = 4*((b12)**2*(A1-A2)**2)/((b12)**2 + (A1-A2)**2)**2
frequencyTCL2 = sp.Rational(1,4)*sp.sqrt(A1**2-2*A1*A2+A2**2+(b12)**2)

modulationDepthTCL2Lambda = sp.lambdify([A1,A2,b12],modulationDepthTCL2,'numpy')
frequencyTCL2Lambda = sp.lambdify([A1,A2,b12],frequencyTCL2,'numpy')

I1o2I1o2WtTCL2Lambda = sp.lambdify([A1,A2,b12,t],I1o2I1o2WtTCL2,'numpy')
I1o2I1o2WtTCL4Lambda = sp.lambdify([A1,A2,b12,t],I1o2I1o2WtTCL4,'numpy')

def I1o2I1o2coherenceTCL2(imap,timeSpace,Bext,highestN,methyl=False,bath=False,mmap={},distMap={}):
    
    As = {}
    Bs = {}
    Js = {}
    
    modulationDepthFrequencyMap = {}
    
    modulationDepthFrequencyMap['m-m pairs'] = {}
    
    modulationDepthFrequencyMap['m-s pairs'] = {}
    modulationDepthFrequencyMap['m-s pairs']['Modulation Depth'] = {}
    modulationDepthFrequencyMap['m-s pairs']['Frequency'] = {}
    modulationDepthFrequencyMap['m-s pairs']['Distances'] = {}
    
    modulationDepthFrequencyMap['s-s pairs'] = {}
    
    ssMods = []
    ssFreqs = []
    ssDists = []
    sseDists = []
    for elem in list(imap.keys()): 
        
        i = elem[0]
        j = elem[1]
    
        if i == 'e':
            As[j] = imap[(i,j)][2][2]
        elif i != 'e':
            Bs[(i,j)] = imap[(i,j)][2][2]
            key = (int(i[:-1]),int(j[:-1]))

            if key in list(mmap.keys()) and methyl:
                Js[(i,j)] = mmap[key]
            else:
                Js[(i,j)] = 0

            
            numI = int(re.match(r"(\d+)([a-zA-Z]+)", i)[1])
            numJ = int(re.match(r"(\d+)([a-zA-Z]+)", j)[1])
            
            typeI = re.findall(r'[a-zA-Z]+', i)[0]
            typeJ = re.findall(r'[a-zA-Z]+', j)[0]
            
            if numI <= highestN and numJ > highestN and typeI != 'F':
        
                modulationDepthFrequencyMap['m-s pairs']['Modulation Depth'][i] = []
                modulationDepthFrequencyMap['m-s pairs']['Frequency'][i] = []
                modulationDepthFrequencyMap['m-s pairs']['Distances'][i] = []
    
    argument = np.zeros(len(timeSpace))
    
    for bint in list(Bs.keys()): 
        k = bint[0]
        l = bint[1]
            
        typeK = re.findall(r'[a-zA-Z]+', k)[0]
        typeL = re.findall(r'[a-zA-Z]+', l)[0]
        
        numK = int(re.match(r"(\d+)([a-zA-Z]+)", k)[1])
        numL = int(re.match(r"(\d+)([a-zA-Z]+)", l)[1])

        if bath:
            if typeK == typeL:
            
                if numL <= highestN: # this condition is sufficent to limit to m-m pairs
                
                    argument += I1o2I1o2WtTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k,l)],timeSpace)

                    modulationDepthFrequencyMap['m-m pairs'][(k,l)] = (modulationDepthTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k,l)]),frequencyTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k,l)]))
                    
                elif numK > highestN: # this condition is sufficient to limit to s-s pairs
                    
                    if distMap[(k,l)] < 3.0 and distMap[(k,l)] > 1.0:

                        argument += I1o2I1o2WtTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k,l)],timeSpace)

                        ssMods.append(modulationDepthTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k,l)]))
                        ssFreqs.append(frequencyTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k,l)]))
                        ssDists.append(distMap[(k,l)])
                        sseDists.append(1/2 * (distMap[('e',k)] + distMap[('e',l)]))
                    
                elif numK <= highestN and numL > highestN: # this condition is sufficient to limit to m-s pairs
                    
                    argument += I1o2I1o2WtTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k, l)],timeSpace)

                    modulationDepthFrequencyMap['m-s pairs']['Modulation Depth'][k].append(modulationDepthTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k,l)]))
                    modulationDepthFrequencyMap['m-s pairs']['Frequency'][k].append(frequencyTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k,l)]))
                    modulationDepthFrequencyMap['m-s pairs']['Distances'][k].append(distMap[(k,l)])
        else:
            if typeK == typeL:
            
                if numL <= highestN: # this condition is sufficent to limit to m-m pairs
                
                    argument += I1o2I1o2WtTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k,l)],timeSpace)

                    modulationDepthFrequencyMap['m-m pairs'][(k,l)] = (modulationDepthTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k,l)]),frequencyTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k,l)]))
            
    # code to sort the three s-s arrays in the same order here

    ssIndices = np.argsort(sseDists)

    ssDistsSorted = np.array(ssDists)[ssIndices]
    sseDistsSorted = np.array(sseDists)[ssIndices]
    ssModsSorted = np.array(ssMods)[ssIndices]
    ssFreqsSorted = np.array(ssFreqs)[ssIndices]
    
    modulationDepthFrequencyMap['s-s pairs']['Modulation Depth'] = ssModsSorted
    modulationDepthFrequencyMap['s-s pairs']['Frequency'] = ssFreqsSorted
    modulationDepthFrequencyMap['s-s pairs']['Distances'] = ssDistsSorted
    modulationDepthFrequencyMap['s-s pairs']['E Distances'] = sseDistsSorted
    
    # code to sort the three m-s pair arrays in the same order
                                                                            
    for key in list(modulationDepthFrequencyMap['m-s pairs']['Distances'].keys()):
        msIndicesKey = np.argsort(modulationDepthFrequencyMap['m-s pairs']['Distances'][key])
        modulationDepthFrequencyMap['m-s pairs']['Distances'][key] = np.array(modulationDepthFrequencyMap['m-s pairs']['Distances'][key])[msIndicesKey]
        modulationDepthFrequencyMap['m-s pairs']['Modulation Depth'][key] = np.array(modulationDepthFrequencyMap['m-s pairs']['Modulation Depth'][key])[msIndicesKey]
        modulationDepthFrequencyMap['m-s pairs']['Frequency'][key] = np.array(modulationDepthFrequencyMap['m-s pairs']['Frequency'][key])[msIndicesKey]

    return 1/2*np.exp(-argument),modulationDepthFrequencyMap

def I1o2I1o2coherenceAPPA(imap,timeSpace,Bext,highestN,methyl=False,bath=False,mmap={},distMap={}):
    
    As = {}
    Bs = {}
    Js = {}
    
    modulationDepthFrequencyMap = {}
    
    modulationDepthFrequencyMap['m-m pairs'] = {}
    
    modulationDepthFrequencyMap['m-s pairs'] = {}
    modulationDepthFrequencyMap['m-s pairs']['Modulation Depth'] = {}
    modulationDepthFrequencyMap['m-s pairs']['Frequency'] = {}
    modulationDepthFrequencyMap['m-s pairs']['Distances'] = {}
    
    modulationDepthFrequencyMap['s-s pairs'] = {}
    
    ssMods = []
    ssFreqs = []
    ssDists = []
    sseDists = []
    for elem in list(imap.keys()): 
        
        i = elem[0]
        j = elem[1]
    
        if i == 'e':
            As[j] = imap[(i,j)][2][2]
        elif i != 'e':
            Bs[(i,j)] = imap[(i,j)][2][2]
            key = (int(i[:-1]),int(j[:-1]))

            if key in list(mmap.keys()) and methyl:
                Js[(i,j)] = mmap[key]
            else:
                Js[(i,j)] = 0

            
            numI = int(re.match(r"(\d+)([a-zA-Z]+)", i)[1])
            numJ = int(re.match(r"(\d+)([a-zA-Z]+)", j)[1])
            
            typeI = re.findall(r'[a-zA-Z]+', i)[0]
            typeJ = re.findall(r'[a-zA-Z]+', j)[0]
            
            if numI <= highestN and numJ > highestN and typeI != 'F':
        
                modulationDepthFrequencyMap['m-s pairs']['Modulation Depth'][i] = []
                modulationDepthFrequencyMap['m-s pairs']['Frequency'][i] = []
                modulationDepthFrequencyMap['m-s pairs']['Distances'][i] = []
    
    argument = np.ones(len(timeSpace))
    
    for bint in list(Bs.keys()): 
        k = bint[0]
        l = bint[1]
            
        typeK = re.findall(r'[a-zA-Z]+', k)[0]
        typeL = re.findall(r'[a-zA-Z]+', l)[0]
        
        numK = int(re.match(r"(\d+)([a-zA-Z]+)", k)[1])
        numL = int(re.match(r"(\d+)([a-zA-Z]+)", l)[1])

        if bath:
            if typeK == typeL:
            
                if numL <= highestN: # this condition is sufficent to limit to m-m pairs
                
                    argument *= (1/2)-2*I1o2I1o2WtTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k,l)],timeSpace)

                    modulationDepthFrequencyMap['m-m pairs'][(k,l)] = (modulationDepthTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k,l)]),frequencyTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k,l)]))
                    
                elif numK > highestN: # this condition is sufficient to limit to s-s pairs
                    
                    if distMap[(k,l)] < 3.0 and distMap[(k,l)] > 1.0:

                        argument *= (1/2)-2*I1o2I1o2WtTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k,l)],timeSpace)

                        ssMods.append(modulationDepthTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k,l)]))
                        ssFreqs.append(frequencyTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k,l)]))
                        ssDists.append(distMap[(k,l)])
                        sseDists.append(1/2 * (distMap[('e',k)] + distMap[('e',l)]))
                    
                elif numK <= highestN and numL > highestN: # this condition is sufficient to limit to m-s pairs
                    
                    argument *= (1/2)-2*I1o2I1o2WtTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k, l)],timeSpace)

                    modulationDepthFrequencyMap['m-s pairs']['Modulation Depth'][k].append(modulationDepthTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k,l)]))
                    modulationDepthFrequencyMap['m-s pairs']['Frequency'][k].append(frequencyTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k,l)]))
                    modulationDepthFrequencyMap['m-s pairs']['Distances'][k].append(distMap[(k,l)])
        else:
            if typeK == typeL:
            
                if numL <= highestN: # this condition is sufficent to limit to m-m pairs
                
                    argument *= (1/2)-2*I1o2I1o2WtTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k,l)],timeSpace)

                    modulationDepthFrequencyMap['m-m pairs'][(k,l)] = (modulationDepthTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k,l)]),frequencyTCL2LambdaJ(As[k],As[l],Bs[(k,l)],Js[(k,l)]))
            
    # code to sort the three s-s arrays in the same order here

    ssIndices = np.argsort(sseDists)

    ssDistsSorted = np.array(ssDists)[ssIndices]
    sseDistsSorted = np.array(sseDists)[ssIndices]
    ssModsSorted = np.array(ssMods)[ssIndices]
    ssFreqsSorted = np.array(ssFreqs)[ssIndices]
    
    modulationDepthFrequencyMap['s-s pairs']['Modulation Depth'] = ssModsSorted
    modulationDepthFrequencyMap['s-s pairs']['Frequency'] = ssFreqsSorted
    modulationDepthFrequencyMap['s-s pairs']['Distances'] = ssDistsSorted
    modulationDepthFrequencyMap['s-s pairs']['E Distances'] = sseDistsSorted
    
    # code to sort the three m-s pair arrays in the same order
                                                                            
    for key in list(modulationDepthFrequencyMap['m-s pairs']['Distances'].keys()):
        msIndicesKey = np.argsort(modulationDepthFrequencyMap['m-s pairs']['Distances'][key])
        modulationDepthFrequencyMap['m-s pairs']['Distances'][key] = np.array(modulationDepthFrequencyMap['m-s pairs']['Distances'][key])[msIndicesKey]
        modulationDepthFrequencyMap['m-s pairs']['Modulation Depth'][key] = np.array(modulationDepthFrequencyMap['m-s pairs']['Modulation Depth'][key])[msIndicesKey]
        modulationDepthFrequencyMap['m-s pairs']['Frequency'][key] = np.array(modulationDepthFrequencyMap['m-s pairs']['Frequency'][key])[msIndicesKey]

    return argument,modulationDepthFrequencyMap

#
#
### Dynamics Functions
#
#

def getCoh(filepath,soi,times,methyl=False,root='',xyzpath=''):
    
    cPos = getCenterSpinDens(filepath) # cPos was originally the set to the origin.
    Bext = np.array([0,0,10000])
    
    xyz,hf = readOrcaV2(filepath,soi)
    
    dim,spins,spin_types = getSpinDictV2(cPos,xyz)
    imap = imapOrcaV2(hf,spins)
    
    highestN = int(re.findall(r'\d+', list(spins.keys())[-1])[0])

    if methyl:
        mmap = methylMap(xyzpath,root)
        return I1o2I1o2coherenceTCL2(imap,times*2*np.pi,Bext,highestN,True,False,mmap)

    return I1o2I1o2coherenceTCL2(imap,times*2*np.pi,Bext,highestN)

def getCohBath(filepath,soi,times,boxLength,bathType,concentration,methyl=False,root='',xyzpath=''):
    
    cPos = getCenterSpinDens(filepath) # necessary to accurately represent the electron - random bath interactions
    
    Bext = np.array([0,0,10000])
    
    xyz,hf = readOrcaV2(filepath,soi)
    molecularAtoms = list(readOrcaV2(filepath,['H','F','C','O','N'])[0].values())

    dim,spins,spin_types = getSpinDictV2(cPos,xyz)
    imap = imapOrcaV2(hf,spins)

    imapB,highestN,distMap = addRandomBath(boxLength,bathType,concentration,spins,imap,molecularAtoms) 

    if methyl:
        mmap = methylMap(xyzpath,root)
        return I1o2I1o2coherenceTCL2(imapB,times*2*np.pi,Bext,highestN,methyl,True,mmap,distMap)

    return I1o2I1o2coherenceTCL2(imapB,times*2*np.pi,Bext,highestN,methyl,True,{},distMap)

#
#
### Electronic Structure Functions
#
#

def getCenterSpinDens(filepath):
    
    spinDens = {}
    
    # if HIRSHFELD analysis is present, extract spin density from that
    # and then skip the Mulliken analysis
    # if HIRSHFELD not present, do the Mulliken analysis as is
    
    with open(filepath,'r') as f:
        
        hirshCount = 0
        
        HIRSHFELD = False
        MULLIKEN = False
        for line in f:
            if 'HIRSHFELD ANALYSIS' in line:
                HIRSHFELD = True
                hirshCount += 1
            elif 'TIMINGS' in line:
                HIRSHFELD = False     
            elif HIRSHFELD:
                if len(line.split()) == 4:
                    spinDens[f'{line.split()[0]}{line.split()[1]}'] = float(line.split()[-1])
            elif 'MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS' in line and hirshCount == 0:
                MULLIKEN = True
            elif 'MULLIKEN REDUCED ORBITAL CHARGES AND SPIN POPULATIONS' in line:
                MULLIKEN = False
            elif MULLIKEN:
                atomNames = (line.split(':'))[0].split()
                if len(atomNames) == 2:
                    spinNumbers = (line.split(':'))[1].split()
                    spinDens[f'{atomNames[0]}{atomNames[1]}'] = float(spinNumbers[-1])
                
    maxSpinDens = max(spinDens, key=spinDens.get)
    
    xyzCount = 0
    with open(filepath,'r') as f:
        XYZ = False
        for line in f:
            if 'CARTESIAN COORDINATES (ANGSTROEM)' in line:
                XYZ = True
                continue
            elif 'CARTESIAN COORDINATES (A.U.)' in line:
                XYZ = False
                continue
            elif XYZ:
                if len(line.split()) == 4: 
                    keyCheck = f'{xyzCount}{line.split()[0]}'
                    if keyCheck == maxSpinDens:
                        posSpinDens = np.array([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])
                    xyzCount += 1
                
    #print(f'max spin density at atom {maxSpinDens}, with position {posSpinDens}')
    return posSpinDens

def readOrcaV2(filepath,soi,efg=False):
    
    with open(filepath,'r') as f:
        for line in f:
            if 'Program Version' in line:
                verscheck = line.split()
                if verscheck[2] == '5.0.4':
                    hfflag = 'Raw HFC matrix'
                    efgflag = 'Raw EFG matrix'
                elif verscheck[2] == '5.0.3':
                    hfflag = 'Raw HFC matrix'
                    efgflag = 'Raw EFG matrix'
                elif verscheck[2] == '6.0.1':
                    hfflag = 'Total HFC matrix'  
                    efgflag = 'Raw EFG matrix'
                elif verscheck[2] == '6.1.1':
                    hfflag = 'Total HFC matrix'
                    efgflag = 'Raw EFG matrix'
                    
    
    xyzdict = {}
    
    nuccount = 0
    
    # getting coordinates
    
    with open(filepath,'r') as f:
        XYZ = False
        for line in f:
            if 'CARTESIAN COORDINATES (ANGSTROEM)' in line:
                XYZ = True
                continue
            elif 'CARTESIAN COORDINATES (A.U.)' in line:
                XYZ = False
                continue
            elif XYZ:
                if len(line.split()) == 4: 
                    if line.split()[0] in soi:
                        xyzdict[f'{nuccount}{line.split()[0]}'] = np.array([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])
                    nuccount += 1
    
    # getting hyperfines
    
    hfdict = {}    
    efgdict = {}
    
    with open(filepath,'r') as f:
        HF = False
        EFG = False
        for line in f:
            if 'Nucleus ' in line and 'Finite Nucleus Model' not in line:
                
                key = line.split(':')[0].split()[1]
                nuc = re.findall(r'[a-zA-Z]+', key)[0]
                if efg:
                    qval = float((next(f)).split()[-2])*1000
                if nuc not in soi:
                    HF = False
                    EFG = False
                    continue
                elif nuc in soi:
                    HF = True
                    if efg:
                        EFG = True
                    continue
            elif 'A(FC)' in line:
                HF = False
                continue
            elif 'V(Tot)' in line:
                if efg:
                    EFG = False
                continue
            elif HF:
                if hfflag in line:
                    hfarr = np.zeros([3,3])
                    dim2count = 0
                    for i in range(4):
                        next_ = next(f).split()
                        if len(next_) == 3:
                            for j,elem in enumerate(next_): 
                                hfarr[dim2count,j] = 1000*float(elem)
                            dim2count += 1
                    hfdict[('0',key)] = hfarr
            elif EFG:
                if efgflag in line:
                    efgarr = np.zeros([3,3])
                    dim2count = 0
                    for i in range(4):
                        next_ = next(f).split()
                        if len(next_) == 3:
                            for j,elem in enumerate(next_): 
                                efgarr[dim2count,j] = float(elem) # no factor of 1000 here since the other function does that part
                            dim2count += 1
                    efgdict[key] = EFGtoQ(efgarr,qval) # transform raw tensor into quadrupolar tensor in kHz, qval should be in millibarn
                    
    
    if not efg:
        return xyzdict,hfdict
    elif efg:
        return xyzdict,hfdict,efgdict

def getSpinDictV2(cPos,xyzs):

    sPos = [cPos]
    spin_types = []
    
    spins = {'e':[cPos,const.gyro_dict['e']]}    
    for spin,key in enumerate(list(xyzs.keys())):
        spins[key] = [xyzs[key],const.gyro_dict[seperate_string_number(key)[-1]]]
        spin_types.append(seperate_string_number(key)[-1])
        
    spin_types.insert(0,'e')
    
    dim = 1
    for type_ in spin_types:
        
        dim *= 2*const.spin_type_dict[type_] + 1
        
    return int(dim), spins, spin_types

def imapOrcaV2(hfobj,spins):
    
    mu0 = np.pi*4e-7 # N A-2
    h = 6.626e-34 # J s
    
    intes = list(spins.keys())
    interactions = combinations(intes,2)
    
    imap = {}
    
    for inte_ in interactions:
        
        i = inte_[0]
        j = inte_[1]
        
        if i == 'e':
            
            imap[(i,j)] = hfobj[('0',j)]
            
        elif i != 'e':
            
            dist = np.sqrt((spins[i][0][0]-spins[j][0][0])**2 + (spins[i][0][1]-spins[j][0][1])**2 + (spins[i][0][2]-spins[j][0][2])**2)
            r12 = spins[i][0] - spins[j][0]
        
            prefactor = ((mu0*spins[i][1]*spins[j][1]*h**2) / (16*np.pi**3*dist**3))*(1000/(h*(1e-10)**3*(1e-4)**2))
            
            aTens = np.array([[(1 - (3*r12[0]**2/dist**2)),-3*r12[0]*r12[1]/dist**2,-3*r12[0]*r12[2]/dist**2],
                    [-3*r12[1]*r12[0]/dist**2,(1 - (3*r12[1]**2/dist**2)),-3*r12[1]*r12[2]/dist**2],
                    [-3*r12[2]*r12[0]/dist**2,-3*r12[2]*r12[1]/dist**2,(1 - (3*r12[2]**2/dist**2))]])
            
            imap[(i,j)] = prefactor*aTens 
            
    return imap


