import os

ROOT = os.path.join(os.path.dirname(__file__), '../')
example_dir = os.path.join(ROOT, 'examples')

for script in os.listdir(example_dir):
    if script.endswith('py'):
        os.system(f"python {example_dir}/{script}")
