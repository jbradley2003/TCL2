import os
import subprocess

ROOT = os.getcwd()

if __name__ == "__main__":
    for d in os.scandir(ROOT):
        if os.path.isdir(d.path) and 'Me' in d.name:
            batch = [f for f in os.listdir(f'{d.path}/scan') if f.endswith('.batch')]
            print(batch)
            for b in batch:
                subprocess.run(f'sbatch run -p msismall {b}',shell=True)