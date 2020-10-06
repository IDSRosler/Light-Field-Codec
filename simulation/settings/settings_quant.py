
import os

# Set required environment variables required for settings.py
os.environ.setdefault('DATASET_DIR', '/home/bosin/lfcodec/datasets')
os.environ.setdefault('RESULTS_DIR', '/home/bosin/lfcodec/results/SIM_QUANT')
os.environ.setdefault('BINARY', 'LF_Codec_QUANT')

from settings import *
