#!/usr/bin/env python

import os
import stat
import sys

# Based on http://stackoverflow.com/a/616686/92396
class Tee(object):
    def __init__(self, name, mode, console=True):
        self.file = open(name, mode)
        self.stdout = sys.stdout
        self.stderr = sys.stderr
        self.console = console
        sys.stdout = self
        sys.stderr = self

    def free(self):
        sys.stdout = self.stdout
        sys.stderr = self.stderr
        self.file.close()

    def write(self, data):
        self.file.write(data)
        if self.console:
            self.stdout.write(data)
        self.flush()

    def flush(self):
        self.file.flush()
        self.stdout.flush()


class spades_error:
    def __init__(self, code, err_str=""):
        self.code = code
        self.err_str = err_str

    def what(self):
        return "Code: " + str(self.code) + " Description: " + self.err_str


def error(err_str, prefix="== Error == ", code=1):
    raise spades_error(code, "\n\n" + prefix + " " + err_str + "\n\n")

#TODO: error log -> log
#TODO: os.sytem gives error -> stop

def sys_call(cmd, cwd=None):
    import shlex
    import time
    import sys
    import subprocess

    cmd_list = shlex.split(cmd)

    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=cwd)

    while not proc.poll():
        if proc.returncode is not None:
            break
        sys.stdout.write(proc.stdout.readline())
        time.sleep(0)

    for line in proc.stdout.readlines():
        print line,
    proc.communicate()

    if proc.returncode:
        error('system call for: "%s" finished abnormally, err code: %d' % (cmd, proc.returncode))


def question_with_timer(question, seconds, default='y'):
    import time
    import curses

    # The timer class    
    class Timer:
        def __init__(self):
            self.target = time.time() + seconds

        def get_left(self):
            return int(self.target - time.time())

    answer = default

    if default == 'n':
        default_str = "[y/N]"
    else:
        default_str = "[Y/n]"
    question[-1] += " " + default_str

    print("\n")
    for line in question:
        print(line)

    t = Timer()
    stdscr = curses.initscr()
    stdscr.nodelay(True)
    curses.noecho()

    try:
        try:
            while True:
                for id, line in enumerate(question):
                    stdscr.addstr(id, 0, line)
                left = t.get_left()
                if left <= 0:
                    break
                stdscr.addstr(len(question), 0, "Seconds left: %02d " % left)
                c = stdscr.getch()
                if c == ord('y'):
                    answer = 'y'
                    break
                elif c == ord('n'):
                    answer = 'n'
                    break
                elif c == ord('\n'):
                    break
        finally:
            # Final operations start here
            stdscr.keypad(0)
            try:
                curses.echo()
                curses.endwin()
            except curses.error, err:
                print("Curses error:", err, "(maybe you are redirecting script's output)")
    except KeyboardInterrupt:
        print("  Exception KeyboardInterrupt was raised: default answer was choosen")

    print("Answer '" + answer + "' was choosen")
    return answer


def save_data_to_file(data, file):
    output = open(file, 'wb')
    output.write(data.read())
    output.close()
    os.chmod(file, stat.S_IWRITE | stat.S_IREAD | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
