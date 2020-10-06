# Argument that are used to create simulation objects.
#
# They are resolved as follows:
#
#   1. Numeric values are inlined.
#   2. List and Tuples are iterated over, creating a
#       separate simulation object for each value.
#   3. Strings are evaluated as previous keys of ARGS.
#       If it fails to find a key, the value is inlined.
#
# Example
# ARGS = {
#     'x': 10,                        # Evaluates to "x 10"
#     'y': [10, 20, 30],              # Create one simulation for
#                                     # each value of y
#     'logfile': '~/mylog.txt'        # As there is no previous entry with
#                                     # the key `logfile`, the string
#                                     # '~/mylog.txt' is inlined.
#                                     # It evaluates to "logfile ~/mylog.txt"
#     'other_x': 'x',                 # Because `x` is a previous key,
#                                     # 'other_x' will have the same value
#                                     # of x.
#                                     # It evaluates to "other_x 10"
# }
ARGS = dict()


# Sequencial arguments to be passed to the executable.
# To set keyword arguments, use `ARGS`
# Example:
# SEQ_ARGS = [ '-one', '-two' ]
SEQ_ARGS = list()

# Executable
EXECUTABLE = None

# Redirect standard input
STDIN = None

# Redirect standard output
STDOUT = None

# Redirect standard error
STDERR = None


# Number of parallel instaces used on the simulation.
# Internally, the implementation uses subprocesses to
# spawn instances. Avoid using more than the number of
# cpus available
PARALLEL_INSTANCES = 1
