import os
import sys

spades_home = os.path.abspath(sys.path[0])
spades_version = ''

def init():
    global spades_home
    global spades_version
    if spades_home == "/usr/bin":
        spades_home = "/usr/share/spades"

    sys.path.append(os.path.join(spades_home, "src/tools/spades_pipeline/"))
    sys.path.append(os.path.join(spades_home, "src/tools/quality/"))
    sys.path.append(os.path.join(spades_home, "src/tools/quality/libs"))

    spades_version =  open(os.path.join(spades_home, 'VERSION'), 'r').readline()
