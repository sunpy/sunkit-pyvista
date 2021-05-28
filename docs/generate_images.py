import os

ROOT = os.path.join(os.path.dirname(__file__), '../')

for script in os.listdir(os.path.join(ROOT, 'examples')):
    if script.endswith('py'):
        os.system(f"python ../examples/{script}")
