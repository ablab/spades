import sys
import re

names = set({"MapSequenceBWA","Substrs", "SHWDistance", "QueueUpdate", "NWDistance", "GetReadAlignment", "BWA", "BestScoredPath", "BestScoredPath_D", "GetWeightedColors", "ProcessCluster"}  )
timetotal = {}
for n in names:
    timetotal[n] = 0
logfile = sys.argv[1]
k = 0
with open(logfile, "r") as log:
    for ln in log:
        for n in names:
            for m in re.finditer("TIME." + n + "=", ln):
                start = m.start()
                #print n, ln[start + len("TIME." + n + "="):].split(" ")[0], float(ln[start + len("TIME." + n + "="):].split(" ")[0])
                timetotal[n] += float(ln[start + len("TIME." + n + "="):].split(" ")[0])
                k += 1
        if k == 10000:
            for n in sorted(names):
                print n, timetotal[n]
            k = 0
            print ""


# names = set({"Substrs", "SHWDistance", "QueueUpdate", "NWDistance", "GetReadAlignment", "BWA", "BestScoredPath", "BestScoredPath_D", "GetWeightedColors", "ProcessCluster"}  )
# timetotal = {}
# for n in names:
#     timetotal[n] = 0
# logfile = sys.argv[1]

# with open(logfile, "r") as log:
#     for ln in log:
#         for n in names:
#             for m in re.finditer("TIME." + n + "=", ln):
#                 start = m.start()
#                 #print n, ln[start + len("TIME." + n + "="):].split(" ")[0], float(ln[start + len("TIME." + n + "="):].split(" ")[0])
#                 timetotal[n] += float(ln[start + len("TIME." + n + "="):].split(" ")[0])

# for n in sorted(names):
#     print n, timetotal[n], timetotal[n]/timetotal["GetReadAlignment"]
#
# print "-----------------------"
# for n in sorted(names):
#     if n == "BestScoredPath_D" or n == "BWA" or n == "BestScoredPath" or n == "GetWeightedColors":
#         print n, timetotal[n], timetotal[n]/(timetotal["GetReadAlignment"] - timetotal["BestScoredPath_D"])
#     else:
#         print n, timetotal[n] - timetotal["BestScoredPath_D"], (timetotal[n] - timetotal["BestScoredPath_D"])/ (timetotal["GetReadAlignment"] - timetotal["BestScoredPath_D"])