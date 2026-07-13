from parseXYZ import *
import subprocess
import os

KEYS = '! B3LYP def2-TZVP def2/J RIJCOSX TightSCF'

def createScan(root, ref_filepaths, filename, charge, mult, group, mpi=False, keywords=KEYS):
  if mpi:
    with open(ref_filepaths[1], 'r') as file:
      lines = file.readlines()
  else:
    with open(ref_filepaths[0], 'r') as file:
      lines = file.readlines()

  params = [
    f'CHARGE = {charge}\n',
    f'MULT = {mult}\n',
    f'AXIS_ATOM1 = {group[0]+1}\n',
    f'AXIS_ATOM2 = {group[1]+1}\n',
    f'ROTOR_H = [{group[2]+1},{group[3]+1},{group[4]+1}]\n',
    f"ORCA_KEYWORDS = '{KEYS}'\n",
]
  lines[27:33] = params
  
  with open(f'{root}/{filename}.py','w') as target:
    target.writelines(lines)
    print(f'{filename}.py has been created.')
    
def createBatch(ref_filepath, dir, input, ntasks, mem):
  with open(ref_filepath,'r') as file:
    lines = file.readlines()

  params = [
    f'#SBATCH --job-name={dir}{input}\n',
    f'#SBATCH --output={input}.log\n',
    f'#SBATCH --ntasks={ntasks}\n',
    f'#SBATCH --mem={mem}\n'
  ]
  lines[1:5] = params
  lines[19] = f'cp {input}.inp $SCRA\n'
  lines[25] = f'/common/software/install/manual/orca/6.1.1/orca {input}.inp > {input}.out\n'

  with open(f'{dir}/scan/{input}.batch','w') as target:
    target.writelines(lines)

def molToScans(mol, ref_filepaths, charge, mult, mem, keywords=KEYS):
  group_map, methyl_pairs = findMethyls(mol)
  
  for central_carbon,group in group_map.items():
    subprocess.run(f'mkdir Me{central_carbon}',shell=True)
    createScan(f'{os.getcwd()}/Me{central_carbon}', ref_filepaths, f'Me{central_carbon}', charge, mult, group, mem, keywords)

def xyzToScans(filename, ref_filepaths, charge, mult, mem, keywords=KEYS):
  mol,atoms,xyz = xyzToMol(filename)
  molToScans(mol, ref_filepaths, charge, mult, mem, keywords)

