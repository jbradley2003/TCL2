import os
import subprocess
from scans import *

ROOT = os.getcwd()
REF = [f'{ROOT}/res/ref_rot.py',f'{ROOT}/res/ref_rot_mpi.py']
XYZ = 'NiL.xyz'
CHARGE = -1
MULT = 2
NSTEPS = 35
NTASKS = 16
MEM = '5g'

KEY = '! B3LYP def2-SVP AUTOAUX RIJCOSX TightSCF'

if __name__ == "__main__":
    xyzToScans(XYZ,REF,CHARGE,MULT,True,KEY)
    for d in os.scandir(os.getcwd()):
        if os.path.isdir(d.path) and 'Me' in d.name:
            py = [f for f in os.listdir(d.path) if f.endswith('.py')]
            command = f'python3 {d.name}/{py[0]} {XYZ} --nsteps {NSTEPS} --ntasks {NTASKS} --mem {MEM}'
            subprocess.run(command, shell=True)
    for d in os.scandir(os.getcwd()):
        if os.path.isdir(d.path) and 'Me' in d.name:
            for i in os.scandir(f'{d.path}/scan'):
                filename = i.name[:-4]
                createBatch(f'{ROOT}/SP.batch',d.name,filename,NTASKS,MEM)
    

