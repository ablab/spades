
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

    def fileno(self):
        if self.stream is not None:
            return self.stream.fileno()

        return self.file.fileno()

def error(err_str, prefix="== Error == ", code=1):
    print("\n\n" + prefix + " " + err_str + "\n\n")
    exit(code)


def sys_call(output_to_console, log_file, cmd):
    import os
    import subprocess
    import sys

    if output_to_console:   cmd += " | tee -a " + log_file
    else:                   cmd += " >> " + log_file

    code = os.system(cmd)

    if code != 0:
        error("command ".join(cmd_list) + " finished abnormally, err code:" + str(code))
