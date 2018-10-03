import os
import subprocess
import sys
import datetime


if len(sys.argv) < 4:
    print "Usage: run_build.py <make dir> <feature> <exe-prefix>"
    exit(-1)
makedir = sys.argv[1]
feature_name = sys.argv[2]
prefix = sys.argv[3]
makecmd = "make -C build/release/projects/segal/ -j4 longreads_aligner"
if not os.path.exists("/home/tdvorkina/results/" + feature_name + "/bin"):
    os.makedirs("/home/tdvorkina/results/" + feature_name + "/bin")
now = str(datetime.datetime.now())
now = now.split(".")[0]
now = now.replace(" ", "_").replace(":", "-")
buildlog = "/home/tdvorkina/results/" + feature_name + "/bin/" + prefix + "_" + now + ".log"
savetolog = makedir + "\n" + makecmd + "\n"
curdir = os.getcwd()

print "Move to " + makedir
os.chdir(makedir)

print "Run git status"
process = subprocess.Popen("git status", shell=True, stdout=subprocess.PIPE, stderr = subprocess.PIPE)
process.wait()
output, error = process.communicate()
if process.returncode:
    print error
    exit(-1)
savetolog += output + "\n"

print "Run git diff"
process = subprocess.Popen("git diff", shell=True, stdout=subprocess.PIPE, stderr = subprocess.PIPE)
process.wait()
output, error = process.communicate()
if process.returncode:
    print error
    exit(-1)
savetolog += output + "\n"

print "Run make"
process = subprocess.Popen(makecmd, shell=True, stdout=subprocess.PIPE, stderr = subprocess.PIPE)
output, error = process.communicate()
if process.returncode:
    print error
    exit(-1)

print "Move exe to /home/tdvorkina/results/" + feature_name + "/bin/"
process = subprocess.Popen("mv ./build/release/bin/longreads_aligner  /home/tdvorkina/results/"  + feature_name + "/bin/" + prefix + "_" + now, shell=True, stdout=subprocess.PIPE, stderr = subprocess.PIPE)
output, error = process.communicate()
if process.returncode:
    print error
    exit(-1)

print "Move to " + curdir
os.chdir(curdir)

fout = open(buildlog, "w")
fout.write(savetolog)
fout.close()
