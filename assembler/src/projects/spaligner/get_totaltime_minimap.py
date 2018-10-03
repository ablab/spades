import sys
import re

# names = set({"BWA", "Minimap"}  )
# timetotal = {}
# for n in names:
#     timetotal[n] = 0
# logfile = sys.argv[1]
# with open(logfile, "r") as log:
#     for ln in log:
#         for n in names:
#             for m in re.finditer(n + " time=", ln):
#                 start = m.start()
#                 #print n, ln[start + len("TIME." + n + "="):].split(" ")[0], float(ln[start + len("TIME." + n + "="):].split(" ")[0])
#                 timetotal[n] += float(ln[start + len(n + " time="):].split(" ")[0])

# for n in names:
#     print n, timetotal[n]

times = {"1m": 60, "5m": 300, "10m": 600}
times_cnt = {"1m": 0, "5m": 0, "10m": 0}
logfile = sys.argv[1]
timetotal = 0

with open(logfile, "r") as log:
    for ln in log:
        for m in re.finditer("read_time=", ln):
            start = m.start()
            #print n, ln[start + len("TIME." + n + "="):].split(" ")[0], float(ln[start + len("TIME." + n + "="):].split(" ")[0])
            t = float(ln[start + len("read_time="):].split(" ")[0])
            timetotal += t
            one = False
            for c in times.keys():
            	if times[c] < t:
            		times_cnt[c] += t
            		if not one:
            			print ln
            			one = True

print "Total=", timetotal
for c in times_cnt.keys():
	print c, times_cnt[c]

               
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