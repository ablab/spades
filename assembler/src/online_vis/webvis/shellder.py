import exceptions

class Shellder:
    #Initialized for sending commands to the vis, reading output, or both
    def __init__(self, pipe_in = None, pipe_out = None, end_out = "[end]\n"):
        if pipe_out is not None:
            self.out = open(pipe_out, "r+")
            self.end_out = end_out
        if pipe_in is not None:
            self.pin = open(pipe_in, "w+")

    def get_line(self):
        if self.out is None:
            raise(IOError("output pipe is not provided!"))
        str = self.out.readline()
        print ("Received " + str)
        if str != self.end_out:
            return str
        else:
            return None
        
    #Reads the whole output and returns as a list of strings
    def get_output(self):
        res = []
        str = self.get_line()
        while str is not None:
            res.append(str)
            str = self.get_line()
        return res

    #Sends a command, like "help"
    def send(self, command):
        if self.pin is None:
            raise(IOError("input pipe is not provided!"))
        print("Sending", command, "...")
        self.pin.write(command + "\n")
        self.pin.flush()
        print("Done")

