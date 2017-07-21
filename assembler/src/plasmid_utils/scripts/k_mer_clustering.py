import sys, os 
import networkx as nx

cutoff = 0.05 

infile = sys.argv[1]

res = os.system ("/Nancy/mrayko/Libs/mash-Linux64-v2.0/mash sketch -i " + infile)
if res != 0:
     print ("mash sketch failed")
     exit(1) 

res = os.system ("/Nancy/mrayko/Libs/mash-Linux64-v2.0/mash dist " + infile + ".msh " + infile + ".msh   | awk '$3<0.1 {print $0}' > " + infile + ".tab")
if res != 0:
     print ("mash dist failed")
     exit(1) 


with open(infile + ".tab", "r") as f:
      pl_table=f.readlines()
 
pl_table = [i.split() for i in pl_table] 

G = nx.Graph()

for i in pl_table:
   if float(i[2]) <= cutoff:
        G.add_edge(i[0], i[1])

for i in nx.connected_components(G):
   print (list(i)[0] + "\t" +  str(len(list(i))))




