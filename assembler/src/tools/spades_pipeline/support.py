
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

def error(err_str, prefix="== Error == ", code=1):
    print("\n\n" + prefix + " " + err_str + "\n\n")
    exit(code)


#TODO: error log -> log
#TODO: os.sytem gives error -> stop

def sys_call(cmd):

    import shlex
    import time
    import sys
    import subprocess

    cmd_list = shlex.split(cmd)
    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    while not proc.poll():
        sys.stdout.write(proc.stdout.readline())
        time.sleep(0)

        if proc.returncode is not None:
            break

    sys.stdout.writelines(proc.stdout.readlines())
    proc.communicate()

    if proc.returncode != 0:
        error("system call for: \"" + cmd + "\" finished abnormally, err code:" + str(proc.returncode))
