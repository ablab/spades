#!/usr/bin/env python

# workaround on stdout & stderr redirecting
class redirected_stream:
    def __init__(self, file, stream=None):
        self.stream = stream
        self.file   = file

    def write(self, data):

        if self.stream is not None:
            self.stream.write(data)
            self.stream.flush()

        self.file.write(data)
        self.file.flush()

    def writelines(self, data):

        if self.stream is not None:
            self.stream.writelines(data)
            self.stream.flush()

        self.file.writelines(data)
        self.file.flush()


    def fileno(self):
        if self.stream is not None:
            return self.stream.fileno()

        return self.file.fileno()

class spades_error:
    def __init__(self, code, err_str = ""):
        self.code    = code
        self.err_str = err_str

def error(err_str, prefix="== Error == ", code=1):
    raise spades_error(code, "\n\n" + prefix + " " + err_str + "\n\n")

#TODO: error log -> log
#TODO: os.sytem gives error -> stop

def sys_call(cmd, cwd = None):

    import shlex
    import time
    import sys
    import subprocess

    cmd_list = shlex.split(cmd)

    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd = cwd)

    while not proc.poll():
        sys.stdout.write(proc.stdout.readline())
        time.sleep(0)

        if proc.returncode is not None:
            break

    sys.stdout.writelines(proc.stdout.readlines())
    proc.communicate()

    if proc.returncode != 0:
        error("system call for: \"" + cmd + "\" finished abnormally, err code:" + str(proc.returncode))

def sys_call_output(cmd):
    import subprocess
    return subprocess.check_output(cmd, shell=True)

def copy(source, dest):
    sys_call("cp " + source + " " + dest)

def question_with_timer(question, seconds):
    import time
    import curses

    # The timer class    
    class Timer():
        def __init__(self):
            self.target = time.time() + seconds
        def get_left(self):
            return int(self.target-time.time())

    t = Timer()
    stdscr = curses.initscr()
    stdscr.nodelay(True)
    curses.noecho()
    answer = 'y'
    while True:
        for id, line in enumerate(question):
            stdscr.addstr(id, 0, line)
        left = t.get_left()
        if left <= 0:
            break
        stdscr.addstr(len(question), 0, 'Seconds left: %s ' % str(left).zfill(2) + ' [default is y]')
        c = stdscr.getch()
        if c == ord('y') :
            break
        elif c == ord('n') :
            answer = 'n'
            break
    # Final operations start here
    stdscr.keypad(0)
    curses.echo()
    curses.endwin()
    return answer
