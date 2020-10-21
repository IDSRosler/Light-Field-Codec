import os
import time
import random
import logging
import subprocess
from glob import glob
from itertools import count
from multiprocessing import Process, Queue, Manager
from . import metrics


class Simulation:
    _count = count()
    instances = Queue()
    instances_status = {}
    status = None

    def __init__(self, executable, args, stdin=None, stdout=None, stderr=None, **options):
        self.executable = executable
        self.args = args
        self.stdin = stdin
        self.stdout = stdout
        self.stderr = stderr

        
        self.id = next(Simulation._count)

        Simulation.instances.put(self)
        Simulation.instances_status.update({
            self.id: {
                'running': False,
                'done': False,
                'time_start': None,
                'time_end': None,
                'worker_pid': None,
                'instance': self,
            }
        })
        
        for key, val in options.items():
            setattr(self, key, val)

    def parse(self):
        line = f'{self.executable}'
        if self.stdin:
            line += f' <{self.stdin}'
        if self.stdout:
            line += f' >{self.stdout}'
        if self.stderr:
            line += f' 2>{self.stderr}'
        if self.args:
            line += ' ' + ' '.join(f'{k} {v}' for k, v in self.args)
        return line

    @classmethod
    def dispatch_all(cls, parallel_instances=1):
        cls.process_lst = []
        cls.manager = Manager()
        Simulation.status = cls.manager.dict()
        Simulation.status.update(Simulation.instances_status)

        for cpu in range(parallel_instances):
            p = Process(target=cls._run, args=(Simulation.status, cpu, ))
            cls.process_lst.append(p)
            p.start()

    @classmethod
    def join(cls):
        for p in cls.process_lst:
            p.join()

    @classmethod
    def _run(cls, status, cpu_affinity=None):
        if cpu_affinity:
            os.sched_setaffinity(0, [cpu_affinity])

        while True:
            try:
                instance = cls.instances.get_nowait()
            except:
                return

            args = map(str, sum(instance.args, ()))
            dargs = dict(instance.args)
            process_info = [instance.executable, *args]
            ref_folder = dargs.get('-input')
            rec_folder = dargs.get('-output')
            if not instance.use_tui:
                print(instance.parse())

            if not instance.fake:
                # Create -output path if it does not exist
                os.makedirs(rec_folder, exist_ok=True)
                open_files = {}
                for stream in ['stdin', 'stdout', 'stderr']:
                    file_stream = getattr(instance, stream)
                    if file_stream:
                        os.makedirs(os.path.dirname(
                            file_stream), exist_ok=True)
                        open_files.setdefault(stream, open(file_stream, 'w'))

            val = status[instance.id]
            val.update({
                'running': True,
                'time_start': time.time()
            })
            status[instance.id] = val

            if not instance.fake:
                subprocess.call(process_info, **open_files) 
                

            val = status[instance.id]
            
            status[instance.id] = val
            if not instance.fake:
                for f in open_files.values():
                    f.close()

            
            if instance.metrics:
                try:
                    metrics.calculate(ref_folder, rec_folder, *instance.metrics)
                except:
                    # Log exception. For now, just silently ignore it.
                    pass
            try:
                for files in instance.discard:
                    for filename in glob(os.path.join(rec_folder, files)):
                        os.remove(filename)
            except:
                pass

            val = status[instance.id]
            val.update({
                'running': False,
                'done': True,
                'time_end': time.time(),
            })
            status[instance.id] = val
            


def build_simulations(settings, options):
    simulations = []
    for args in settings.args:
        parsed_values, args = settings.parse_variables(args)
        sim = Simulation(settings.EXECUTABLE,
                         list(args.items()),
                         stdin=parsed_values['STDIN'],
                         stdout=parsed_values['STDOUT'],
                         stderr=parsed_values['STDERR'],
                         **vars(options))
        simulations.append(sim)

        # only runs one simulation
        if options.once:
            break
    return simulations
