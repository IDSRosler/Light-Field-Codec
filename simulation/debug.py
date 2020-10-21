# Module used to override the optins set in settings.py
# This is done to create a controlled environment to run very 
# specific configurations. Environment variables are set previuosly 
# to importing the module.

import os

# Set required environment variables required for settings.py
os.environ.setdefault('DATASET_DIR', '/mnt/c/Users/Cristian/Documents/UFPel/lfcodec/Datasets')
os.environ.setdefault('RESULTS_DIR', '/mnt/c/Users/Cristian/Documents/UFPel/lfcodec/Results')
os.environ.setdefault('BUILD_DIR', 'build')

from settings import *

# Update ARGS 
ARGS.update({
    '-transform': ['DST', 'DCT'],
    '-qp': [15],
})
