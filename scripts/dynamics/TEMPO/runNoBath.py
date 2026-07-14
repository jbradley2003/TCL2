import os
import numpy as np
import pickle

import sys

from functions import *
from constants import *

# Concentrations

avogadrosNumber = 6.022e23

# (1:1) (w:w) --> @ 1g we have (1/0.9167) mL H20 : (1/1.261) mL Gly (Hahn echo is performed @ 20 K)

densityIce = 0.9167 # g/mL at 0 deg C... - this is a placeholder.
molarMassIce = 18.02
numberProtonsIce = 2

concentrationIce = densityIce / molarMassIce * 1000 / 1000 * avogadrosNumber * numberProtonsIce

densityGly = 1.261 # g/mL at 0 deg C... - this is a placeholder.
molarMassGly = 92.09 
numberProtonsGly = 8

concentrationGly = densityGly / molarMassGly * 1000 / 1000 * avogadrosNumber * numberProtonsGly

# Weighted average using the volume fractions of the two components.
concentrationMixture = (concentrationIce * (1 / densityIce) + concentrationGly * (1 / densityGly)) / ((1 / densityIce) + (1 / densityGly))

# densityDCM = 1.33 # g/mL
# molarMassDCM = 84.93 # g/mol
# numberProtonsDCM = 2

# concentrationDCM = densityDCM / molarMassDCM * 1000 / 1000 * avogadrosNumber * numberProtonsDCM

# densityToluene = 0.867 # g/mL
# molarMassToluene = 92.14 # g/mol
# numberProtonsToluene = 8

# concentrationToleune = densityToluene / molarMassToluene * 1000 / 1000 * avogadrosNumber * numberProtonsToluene

# concentrationMixture = (concentrationDCM*9 + concentrationToleune)/10

concentration = concentrationMixture
boxLength = 30
bathType = 'H'

times = np.linspace(0,0.02,120)

numBathInstances = 1000

averageCoherence = np.zeros([len(times)])

# Paths
hfc_path = './../../ORCA/Hyperfine/tempo.out'

# Data Initialization

numBathSpins,numBathPairs,molecularKeys = initializeStatisticsV2(hfc_path,['H','F'],boxLength,bathType,concentration)

msDistances = {}
msFrequencies = {}
msModulationDepths = {}

# initialize m-s statistics objects

for key in molecularKeys:
    msDistances[key] = np.zeros([numBathSpins])
    msFrequencies[key] = np.zeros([numBathSpins])
    msModulationDepths[key] = np.zeros([numBathSpins])
    
# initialize s-s statistics

ssDistances = np.zeros([numBathPairs])
sseDistances = np.zeros([numBathPairs])
ssFrequencies = np.zeros([numBathPairs])
ssModulationDepths = np.zeros([numBathPairs])


for i in range(numBathInstances):
    
    coherence,modDepthFreqMap = getCoh(hfc_path,['H','F'],times)
    
    averageCoherence += coherence
    
    # average over m-s pairs here
    
    for key in list(modDepthFreqMap['m-s pairs']['Distances'].keys()):
        
        msDistances[key][:len(modDepthFreqMap['m-s pairs']['Distances'][key])] += modDepthFreqMap['m-s pairs']['Distances'][key]
        msFrequencies[key][:len(modDepthFreqMap['m-s pairs']['Frequency'][key])] += modDepthFreqMap['m-s pairs']['Frequency'][key]
        msModulationDepths[key][:len(modDepthFreqMap['m-s pairs']['Modulation Depth'][key])] += modDepthFreqMap['m-s pairs']['Modulation Depth'][key]
    
    # average over s-s pairs here
    
    ssDistances[:len(modDepthFreqMap['s-s pairs']['Distances'])] += modDepthFreqMap['s-s pairs']['Distances']
    sseDistances[:len(modDepthFreqMap['s-s pairs']['E Distances'])] += modDepthFreqMap['s-s pairs']['E Distances']
    ssFrequencies[:len(modDepthFreqMap['s-s pairs']['Frequency'])] += modDepthFreqMap['s-s pairs']['Frequency']
    ssModulationDepths[:len(modDepthFreqMap['s-s pairs']['Modulation Depth'])] += modDepthFreqMap['s-s pairs']['Modulation Depth']
    
    print(f'{i+1} out of {numBathInstances}')
    
averageCoherence /= numBathInstances

# divide each m-s value by numBathInstances here

for key in list(msDistances.keys()):
    
    msDistances[key] /= numBathInstances
    msFrequencies[key] /= numBathInstances
    msModulationDepths[key] /= numBathInstances

# divide each s-s value by numBathInstances here

ssDistances /= numBathInstances
sseDistances /= numBathInstances
ssFrequencies /= numBathInstances
ssModulationDepths /= numBathInstances

# creating dictionaries of these objects

msDictionary = {'Distances':msDistances,
               'Frequency':msFrequencies,
               'Modulation Depth':msModulationDepths}

ssDictionary = {'Distances':ssDistances,
                'E Distances': sseDistances,
               'Frequency':ssFrequencies,
               'Modulation Depth':ssModulationDepths}

outputDictionary = {'Concentration of Spin Bath (spins / cm^3)':concentration,
                   'Number of Bath Orientations in Average':numBathInstances,
                   'Box Edge Length':boxLength,
                   'Times (millisecond)':times,
                   'Coherence':averageCoherence,
                   'm-m Mod. Depth and Frequency':modDepthFreqMap['m-m pairs'],
                   'm-s statistics':msDictionary,
                   's-s statistics':ssDictionary}

with open(f'tempo_bL{boxLength}_nB{numBathInstances}_c10_TCL2_filtered_geminal_no_methyl_no_bath.pkl','wb') as f:
    pickle.dump(outputDictionary,f)
