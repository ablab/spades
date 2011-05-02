import sys
import subprocess 
import gzip

def safeOpenFile(fileToOpen, mode, errorMessage):
    try:
        if ".gz" in fileToOpen:
            openedFile = gzip.open(fileToOpen, mode)
        else:
            openedFile = open(fileToOpen, mode)
        return openedFile
    except IOError:
        sys.exit(errorMessage)

def fileLinesCount(fileName):
    p = subprocess.Popen(['wc', '-l', fileName], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode:
        raise IOError(err)
    return int(result.strip().split()[0]) 

