import os, time
import multiprocessing

from threading import Thread
from threading import Lock

# estrura *CMD*; dict {bin: ;argv: ;output: }

class PlataformaExecucaoMultithread():
    """Parametros comados, numThreads é opicional"""

    def __init__(self, comandos, output_file,  numThreads=5):
        self.plaf_folder = output_file
        self.log =  output_file + 'Logs/'
        self.comandos = comandos
        self.numThreads = numThreads
        self.mutex = Lock()

        if not os.path.isdir(self.plaf_folder):
            os.system('mkdir ' + self.plaf_folder)
        
        if not os.path.isdir(self.log):
            os.system('mkdir ' + self.log)

        for d in self.comandos:
            if not os.path.isdir(self.plaf_folder + d['output']):
                os.system('mkdir ' + self.plaf_folder + d['output'])

        self.start()

    def start(self):
        for i in range(self.numThreads):
            t = Thread(target=self._work, args=[i])
            t.start()

    def _work(self, i):
        while len(self.comandos):

            if len(self.comandos):
                cmd = self.comandos.pop(0)

                if os.path.exists(self.plaf_folder + cmd['output']):
                    
                    curr_cmd = ' '.join(['taskset -c ' + str(i+1), cmd['bin'], cmd['argv'], '>', self.log + cmd['output']])
                    
                    print(curr_cmd)
                    
                    os.system(curr_cmd)

                    os.system("touch " + self.plaf_folder + cmd['output'].split('/')[-1])


if __name__ == "__main__":
    print("main")