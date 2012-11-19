import sys, os


genome_file = open(sys.argv[1])

first_line = True

genome = ""
for line in genome_file:
  if first_line:
    first_line = False
    continue

  genome += line.strip()


graph_file = open(sys.argv[2])
out_file = open(sys.argv[3], "w")


s = set()
ind = 0
for line in graph_file:
  
  if line.find("->") == -1:
    continue

  
  pos = int(line.split("pos=")[1].split(" ")[0])
  length = int(line.split("len=")[1].split(" ")[0])
  ch = line.split("ch=")[1].split("'")[1]
  
  edge = line.split("[")[0].strip().split("->")
  from_v = int(edge[0])
  to_v = int(edge[1])
  
  triple = (from_v, to_v, ch)
  
  if triple in s:
    print "Duplicate", triple
  else:
    out_file.write(">contig" + str(ind) + "\n")
    out_file.write(genome[pos : pos + length] + "\n")
    ind += 1
    s.add(triple)
  
  

out_file.close()
  


