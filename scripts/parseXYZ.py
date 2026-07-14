import py3Dmol
import numpy as np
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw
from itertools import combinations
from xyz2mol import *
import io
import os

ROOT = os.getcwd() + '/ORCA/Barrier/'
NIL = f'{ROOT}NiL.xyz'
TEMPO = f'{ROOT}tempo.xyz'
TM = f'{ROOT}tmQM_X1.xyz'

# Phase 1

def readXYZ(filepath):
    coords = {}
    with open(filepath, 'r') as file:
        for line_no, line in enumerate(file, start=1):
            tmp = line.split()
            if len(tmp) == 4:
                coords[line_no-3] = [tmp[0], np.array(tmp[1:], dtype=float)]
    return coords

def dist(pos1, pos2):
    return np.sqrt(np.sum((pos1-pos2)**2))

def distMatrix(coords):
    atoms = list(coords.keys())
    pairs = list(combinations(atoms,2))
    d = np.zeros((len(atoms),len(atoms)))
    for (i,j) in pairs:
        d[i,j] = dist(coords[i][1], coords[j][1])
        d[j,i] = d[i,j]
    return d

# Depends on bounds set for sp^3 C-H and C-C bond lengths.
def findMethyls(filepath,bondrange=[0,1.15,1.3,1.6]):
    coords = readXYZ(filepath)
    d =  distMatrix(coords)
    carbons = [k for k,v in coords.items() if v[0] == 'C']
    adj_list = []
    for c in carbons:
        distList = d[c,:]
        adjC = []
        adjH = []
        for i in range(len(distList)):
            if distList[i] > bondrange[0] and distList[i] <= bondrange[1]:
                adjH.append(f'{i}{coords[i][0]}')
            if distList[i] > bondrange[2] and distList[i] <= bondrange[3]:
                adjC.append(f'{i}{coords[i][0]}')
        if len(adjH) == 3:
            adj_list.append([f'{c}C', adjC, adjH])
    return adj_list

# Phase 2

def xyzToMol(filepath):
    atoms, charge, xyz = read_xyz_file(filepath)
    # May need to change this to toggle Extended Huckel
    mol = xyz2mol(atoms,xyz,charge,True,True,True,False,False)
    return mol, atoms, xyz

def viewMol3D(mol,name='',w=500,h=500):
    mb = Chem.MolToMolBlock(mol)
    view = py3Dmol.view(width=w, height=h)
    view.addModel(mb, "sdf")
    
    view.setStyle({
        'stick': {'radius': 0.15},  # The "sticks"
        'sphere': {'scale': 0.3}    # The "balls" (scaled down from standard VdW)
    })
    view.zoomTo()
    out_html = os.path.join(os.getcwd(), f"html/{name}_molecule_view.html")
    view.write_html(out_html)
    print(f"Saved 3D view to {out_html}")

def drawMol(mol,w=500,h=500,rot_angle=0):
    img = Draw.MolToImage(mol)
    d2d = Draw.MolDraw2DCairo(w, h)
    dopts = d2d.drawOptions()
    dopts.rotate = rot_angle
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    img = Image.open(io.BytesIO(d2d.GetDrawingText()))
    img.show()

def bondsInSubstructMatch(mol,substruct):
    mol.UpdatePropertyCache()
    pattern = Chem.MolFromSmarts(substruct)
    return mol.GetSubstructMatches(pattern)

# Finds methyl groups based on atomic connectivity determined with methods in xyz2mol package.
def findMethyls(mol):
    # [CH3][*] --> Looking for R-CH3 groups
    methyl_bonds = bondsInSubstructMatch(mol,'[CH3][*]')
    central_carbons = set([i for (i,j) in methyl_bonds])
    parent_atoms = set([j for (i,j) in methyl_bonds if mol.GetAtomWithIdx(j).GetSymbol() != 'H'])
    # Stores methyl groups in the following order [Parent atom, C, H, H, H]
    methyls = {}
    for i in list(central_carbons):
        group = []
        hydros = []
        child = mol.GetAtomWithIdx(i)
        for j in list(parent_atoms):
            if (i,j) in methyl_bonds:
                group.append(j)
                group.append(i)
        for bond in child.GetBonds():
            k = bond.GetEndAtomIdx()
            if k != i:
                group.append(k)
                hydros.append(k)
        methyls[i] = group, combinations(hydros,2)    
    return methyls, methyl_bonds

# Useful for pulling individual structures from a large XYZ file.
def createXYZs(filepath,base_name,n=0):
    file_count = 0
    with open(filepath, 'r') as file:
        for line in file:
            tmp = line.split()
            if len(tmp) == 1:
                index = 0
                atom_count = int(tmp[0])
                file_count += 1
                file_name = f'{base_name}{file_count}.xyz'

                if n > 0 and file_count > n:
                    break
                with open(file_name, 'w') as target:
                    target.write("\n" * (atom_count+3))
                with open(file_name, 'r') as target:
                    lines = target.readlines()       
            lines[index] = line
            index += 1 
            if index == atom_count + 2:
                with open(file_name, 'w') as target:
                    target.writelines(lines)

