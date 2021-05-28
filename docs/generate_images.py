import os
for script in os.listdir('../examples'):
    if script.endswith('py'):
        os.system(f"python ../examples/{script}")
