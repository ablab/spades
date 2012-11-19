import sys

if len(sys.argv) != 4:
  print "<input nucmer report> <input mis_log> <output>"
  exit(1)
tmp_file_name = "information_tmp.txt"
CONTIG = "CONTIG"
MISASSEMBL = "Extensive misassembly"
ALIG = "Real Alignment"
out = open(tmp_file_name, "w")
infos = []
for line in open(sys.argv[1]):
  if CONTIG in line:
    is_mis = False
    for info in infos:
      if MISASSEMBL in info:
        is_mis = True
        break
    if is_mis:
      out.write("-----------------------\n")
      important_info = ""
      align = 0
      for info in infos:
        if CONTIG in info: 
          important_info += info.split("_")[1].strip() + "\n"
        elif ALIG in info:
          important_info += str(align) + ":" + info.split("|")[1].strip() + "\n"
          align +=1
        out.write(info)
      out.write(important_info)
    infos = [line]
  else:
    infos.append(line)
out.write("-----------\n")

out.close()
nucm_info = open(tmp_file_name)

infos = []
all_contigs = dict()
id_cont = 0
rc_id = 0
for line in nucm_info:
  if line[0] == '-':
    if id_cont != rc_id:
      all_contigs[id_cont] = infos
     # all_contigs[rc_id] = infos
    infos = []
  elif CONTIG in line:
    ids_contig = line.split('_')
    id_cont =  int(ids_contig[2])
    rc_id = int(ids_contig[3])
  infos.append(line)
log = open(sys.argv[2])
out = open(sys.argv[3], "w")
MIS = "misassemble"
line = log.readline()
while line:
  if MIS in line:
    lst = line.split()
    id1 = int(lst[1])
    if id1 in all_contigs:
      infos = all_contigs[id1]
      for info in infos:
        out.write(info)
      out.write("\n" + line)
      line = log.readline()
      while line and not MIS in line:
        out.write(line)
        line = log.readline()
    else:
      line = log.readline()
  else:
    line = log.readline()
      
out.close()

