#!/usr/bin/python3
import argparse
import sys
import curses
from core.settings import Settings
from core.simulation import Simulation, build_simulations
from core import screen


PY_VERSION = (3, 6)
if not sys.version_info >= PY_VERSION:
    print("Python %d.%d or higher is required." % PY_VERSION)
    print("You are using Python %d.%d.%d." % sys.version_info[:3])
    sys.exit(1)

parser = argparse.ArgumentParser(description='Best simulator in town!')
parser.add_argument('-f', '--fake', action='store_true',
                    help='Print out the commands but does not run')
parser.add_argument('-x', '--use-tui', action='store_true',
                    help='Experimental Text User Interface')
parser.add_argument('-1', '--once', action='store_true',
                    help='Only runs one simulation and stops.')
parser.add_argument('-S', '--settings', default=None,
                    help='Dotted path to the settings module.')
parser.add_argument('-o', '--display-stdout', action='store_true',
                    help='If enabled, display stdout regardless of '\
                        'whichever redirection defined in settings.')
parser.add_argument('-m', '--metrics', action="extend", nargs="+", type=str,
                    help='Metrics to be calculated after a simulation '\
                        'finishes.')
parser.add_argument('-d', '--discard', action="extend", nargs="+", type=str,
                    help='If set, discards the created .')


def tui(stdscr):
    total_jobs = len(Simulation.instances_status)
    s = screen.SimulationScreen(stdscr, total_jobs)
    s.init()
    s.loop(settings.PARALLEL_INSTANCES)

if __name__ == "__main__":
    options = parser.parse_args()   
    settings_module = options.settings or 'settings'
    settings = Settings(settings_module)
    if options.display_stdout:
        setattr(settings, 'STDOUT', None)
    build_simulations(settings, options=options)
    if options.use_tui:
        curses.wrapper(tui)
    else:
        Simulation.dispatch_all(settings.PARALLEL_INSTANCES)
        Simulation.join()
