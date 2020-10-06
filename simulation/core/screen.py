import curses
import time
import sys
import re
from simulation import Simulation

start_time = int(time.time())


class SimulationScreen:
    COLOR_RUNNING = None
    COLOR_WAITING = None
    COLOR_DONE = None
    COLOR_ERROR = None

    def __init__(self, stdscr, total_jobs):
        self.H, self.W = stdscr.getmaxyx()
        self.stdscr = stdscr
        self.total_jobs = total_jobs
        self.title_bar = curses.newwin(1, self.W, 0, 0)
        self.status_bar = curses.newwin(1, self.W, self.H - 1, 0)
        self.info_pad = curses.newpad(self.total_jobs, 22)
        self.job_pad = curses.newpad(self.total_jobs, 1024)

        self.info_offset = [0, 0]
        self.job_offset = [0, 0]
        self.timers = {}
        self.jobs_finished = 0
        self.jobs_instances = []
        self.shrink_paths = False
        self.graceful_exit = False
        self._end_loop = False
        self.status_message = ''

    def recreate_windows(self):
        self.H, self.W = self.stdscr.getmaxyx()
        self.title_bar = curses.newwin(1, self.W, 0, 0)
        self.status_bar = curses.newwin(1, self.W, self.H - 1, 0)
        self.info_pad = curses.newpad(self.total_jobs, 22)
        self.job_pad = curses.newpad(self.total_jobs, 1024)

        for i, (_, job) in enumerate(Simulation.status.items()):
            self.info_pad.addstr(i, 0, '[                   ]')
            self.info_pad.addstr(
                i, 1, 'AWAITING'.center(19), self.COLOR_WAITING)
            parsed_job = job['instance'].parse()
            if self.shrink_paths:
                parsed_job = ' '.join(f'{path[:6]}~{path[-8:]}'
                                      if re.match(r"([><]?[/\.]+[a-zA-Z\./]*[\s]?)", path) else path
                                      for path in parsed_job.split())
            self.job_pad.addstr(i, 0, parsed_job.ljust(1023))

    def elapsed_time(self, timer):
        start = self.timers.setdefault(timer, time.time())
        return time.strftime('%H:%M:%S', time.gmtime(time.time() - start))

    def sec2str(self, sec):
        return time.strftime('%H:%M:%S', time.gmtime(sec))

    def loop(self, parallel_instances):

        Simulation.dispatch_all(parallel_instances)
        self.recreate_windows()
        self.refresh_all()

        next_c = None
        while not self._end_loop:
            if time.time() - self.last_refresh > 0.25:
                self.refresh_all()
            c = next_c
            next_c = self.stdscr.getch()
            if next_c != -1:
                continue

            if c == ord('q'):
                self.graceful_exit = True
                self.set_status('Cleaning queue')
                while not Simulation.instances.empty():
                    Simulation.instances.get_nowait()
                self.set_status('Waiting processes to finish')
            elif c == ord('s'):
                self.shrink_paths = not self.shrink_paths
                self.update_jobpad()
            elif c == curses.KEY_LEFT:
                self.job_offset[1] = max(self.job_offset[1] - 1, 0)
                self.update_jobpad()
            elif c == curses.KEY_RIGHT:
                self.job_offset[1] = min(self.job_offset[1] + 1, 1024)
                self.update_jobpad()

            elif c == curses.KEY_DOWN:
                height = self.total_jobs - self.H - 2
                self.job_offset[0] = min(self.job_offset[0] + 1, height)
                self.info_offset[0] = min(self.info_offset[0] + 1, height)
                self.update_jobpad()
                self.update_infopad()
            
            elif c == curses.KEY_UP:
                self.job_offset[0] = max(self.job_offset[0] - 1, 0)
                self.info_offset[0] = max(self.info_offset[0] - 1, 0)
                self.update_jobpad()
                self.update_infopad()

            elif c == curses.KEY_NPAGE:
                height = self.total_jobs - self.H - 2
                self.job_offset[0] = min(self.job_offset[0] + 20, height)
                self.info_offset[0] = min(self.info_offset[0] + 20, height)
                self.update_jobpad()
                self.update_infopad()
                
            elif c == curses.KEY_PPAGE:
                self.job_offset[0] = max(self.job_offset[0] - 20, 0)
                self.info_offset[0] = max(self.info_offset[0] - 20, 0)
                self.update_jobpad()
                self.update_infopad()

            elif c == curses.KEY_RESIZE:
                self.recreate_windows()

        self.refresh_all()
        self.set_status('Joining processes and exiting.')
        Simulation.join()

    def update_titlebar(self):
        _, width = self.title_bar.getmaxyx()
        title = ' '.join(sys.argv)

        timestamp = self.elapsed_time('titlebar')
        line = f' {title:<{width-12}}{timestamp:>9} '
        self.title_bar.addstr(0, 0, line, curses.A_REVERSE)
        self.title_bar.refresh()

    def update_statusbar(self):
        _, width = self.status_bar.getmaxyx()
        jobs = f' Jobs: {self.jobs_finished}/{self.total_jobs}'
        line = f'{jobs:<{15}}{self.status_message:>{width - 17}} '
        self.status_bar.addstr(0, 0, line, curses.A_REVERSE)
        self.status_bar.refresh()

    def set_status(self, msg):
        self.status_message = msg
        self.status_bar.refresh()

    def update_jobpad(self):
        status = list(Simulation.status.items())
        for i, (_, job) in enumerate(status):
            parsed_job = job['instance'].parse()
            if self.shrink_paths:
                parsed_job = ' '.join(f'{path[:6]}~{path[-8:]}'
                                      if re.match("([><]?[/\.]+[a-zA-Z\./]*[\s]?)", path) else path
                                      for path in parsed_job.split())
            self.job_pad.addstr(i, 0, parsed_job.ljust(1023))
        self.job_pad.refresh(*[*self.job_offset, 2, 23, self.H - 2, self.W-2])

    def update_infopad(self):
        status = list(Simulation.status.items())
        finished = 0
        running = 0
        for i, (_, job) in enumerate(status):
            if job['running'] and not job['done']:
                line = f' {self.elapsed_time(f"timer{i}"):<8}'
                self.info_pad.addstr(i, 1, line.ljust(19))
                self.info_pad.addstr(i, 11, 'RUNNING', self.COLOR_RUNNING)
                running += 1

            elif not job['running'] and job['done']:
                finished += 1
                seconds = job['time_end'] - job['time_start']
                line = f' {self.sec2str(seconds):<8}'
                self.info_pad.addstr(i, 1, line.ljust(19))
                self.info_pad.addstr(i, 11, ' DONE!', self.COLOR_DONE)

            elif self.graceful_exit:
                self.info_pad.addstr(
                    i, 1, 'CANCELED'.center(19), self.COLOR_ERROR)

        self.jobs_finished = finished
        if running == 0:
            self._end_loop = True
        self.info_pad.refresh(*[*self.info_offset, 2, 1, self.H - 2, 22])

    def refresh_all(self):
        self.update_titlebar()
        self.update_statusbar()
        self.update_jobpad()
        self.update_infopad()
        self.last_refresh = time.time()

    def cleanup(self):
        curses.curs_set(1)
        curses.nocbreak()
        self.stdscr.keypad(False)

        curses.endwin()

    def init(self):
        curses.start_color()
        curses.init_pair(1, curses.COLOR_BLUE, curses.COLOR_BLACK)
        curses.init_pair(2, curses.COLOR_YELLOW, curses.COLOR_BLACK)
        curses.init_pair(3, curses.COLOR_GREEN, curses.COLOR_BLACK)
        curses.init_pair(4, curses.COLOR_RED, curses.COLOR_BLACK)
        curses.noecho()
        curses.cbreak()
        curses.curs_set(0)
        self.stdscr.keypad(True)
        self.stdscr.nodelay(True)

        self.COLOR_RUNNING = curses.color_pair(1)
        self.COLOR_WAITING = curses.color_pair(2)
        self.COLOR_DONE = curses.color_pair(3)
        self.COLOR_ERROR = curses.color_pair(4)
