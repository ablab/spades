import sys

if len(sys.argv) < 3:
  print("<input gap file> <output distribution>")
  exit(1)

fin = open(sys.argv[1])
fout = open(sys.argv[2], "w")

line = fin.readline()
line = fin.readline()
count = 0
distr = dict()
while line:
  lst = line.split(" ")
  begin = int(lst[0])
  end = int(lst[1])
  dist = end - begin
  if dist not in distr:
    distr[dist] = 0
  distr[dist] += 1
  count += 1
  line = fin.readline()

fout.write(str(count) + "\n")
for dist, dist_count in distr.items():
  fout.write(str(dist) + " " + str(dist_count) + "\n")
fout.close()
