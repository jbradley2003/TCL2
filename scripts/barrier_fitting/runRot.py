import os
import subprocess

xyz = 'NiL.xyz'
mem = '60g'
ntasks = 16
nsteps = 40

for d in os.scandir(os.getcwd()):
    if os.path.isdir(d.path) and 'Me' in d.name:
         py = [f for f in os.listdir(d.path) if f.endswith('.py')]
         command = f'python3 {d.name}/{py[0]} {xyz} --nsteps {nsteps} --ntasks {ntasks} --mem {mem}'
         subprocess.run(command, shell=True)
