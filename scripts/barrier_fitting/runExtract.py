import os
import subprocess

for d in os.scandir(os.getcwd()):
    if os.path.isdir(d.path) and 'Me' in d.name:
        command = f'python3 {d.name}/scan extract.py --out {d.name}/barrier.csv'
        subprocess.run(command, shell=True)