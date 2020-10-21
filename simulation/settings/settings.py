# Module to configure simulation parameters. 

import os
from os.path import dirname, abspath, join

# Helper method to read from environment variables.
# If not quiet, it raises an exception when variable is not defined.
def from_env(variable, quiet=False):
    value = os.getenv(variable)
    if value is None and not quiet:
        raise KeyError(f'Environment variable \'{variable}\' not defined')
    return value


# Path to parent directory.
BASE_DIR = dirname(dirname(dirname(abspath(__file__))))


# The path to folder containing all datasets
DATASET_DIR = from_env('DATASET_DIR')

# The path to folder to save all results
RESULTS_DIR = from_env('RESULTS_DIR')

# List of datasets to run on this simulation.
DATASETS_TO_RUN = ['Bikes', ]

# Resolved paths for each dataset.
# The executable expects a trailing slash for folders.
DATASETS = [join(DATASET_DIR, dataset, '')
            for dataset in DATASETS_TO_RUN]
RESULTS = [join(RESULTS_DIR, dataset + '_%(SIMULATION_ID)s', '')
           for dataset in DATASETS_TO_RUN]

# List of transformation protocols to be used
TRANSFORMS = ['DCT']

# Quatizations values
QUANTIZATIONS = [1, 3, 5, 7]

# Axis to flip
AXIS_X = 1 << 0
AXIS_Y = 1 << 1
AXIS_U = 1 << 2
AXIS_V = 1 << 3
FLIPS = [
    0,
    # AXIS_U,
    # AXIS_V,
    # AXIS_U | AXIS_V,
    # AXIS_V | AXIS_X,
    # AXIS_U | AXIS_X,
]

# Segments
NO_SEGMENTS = 0
SIDE_8 = 1
SIDE_4 = 2
SIDE_2 = 3
SEGMENTS = [0, 1, 2]

ARGS = {
    '-transform': TRANSFORMS,
    '-qp': QUANTIZATIONS,
    '-qx': '-qp',       # Reference to -qp
    '-qy': '-qp',       # Reference to -qp
    '-qu': '-qp',       # Reference to -qp
    '-qv': '-qp',       # Reference to -qp
    '-lfx': 625,
    '-lfy': 434,
    '-lfu': 13,
    '-lfv': 13,
    '-blx': 15,
    '-bly': 15,
    '-blu': 13,
    '-blv': 13,
    '-flipaxis': FLIPS,
    '-segments': SEGMENTS,
    '-input': DATASETS,
    '-output': RESULTS,  # This parameter uses SIMULATION_ID
}

SEQ_ARGS = ('-lytro', )


GROUP_TOGETHER_ARGS = (
    ('-input', '-output'),
)

BINARY = from_env('BINARY')

# Binary
EXECUTABLE = join(BASE_DIR, 'build', BINARY)

# Variable at module level.
SIMULATION_ID = '%(-qp)s_%(-transform)s_%(-segments)s_%(-flipaxis)s'

STDOUT = '%(-output)s' + '%(SIMULATION_ID)s.txt'

VARIABLES_TO_PARSE = [
    'SIMULATION_ID',
    '-output',
]

# Number of parallel instaces to run the simulation
PARALLEL_INSTANCES = 6

PLUGINS = [
    'pysim.plugins.LogEvents',
    'pysim.plugins.SaveIntermediaryState',
]
