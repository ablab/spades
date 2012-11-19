import sys

if len(sys.argv) != 5:
  print "<input nucmer file> <begin pos> <end pos> <output>"
  exit(1)

CONTIG = "CONTIG"
UNALIGNED_CONTIG = "This contig is unaligned."
REAL_ALIGN = "Real Alignment"
TOTAL_ALIGN = "One align captures most of this contig"
MISASSEMBl = "Exte"
nucmer = open(sys.argv[1], "r")
out = open(sys.argv[4], "w")
begin_pos = int(sys.argv[2])
end_pos = int(sys.argv[3])
infos = []
result_info = []
for line in nucmer:
  if CONTIG in line:
    if True:
      between_poses = False
      first_pos = 100000000000000000000000
      important_info = "-----------------\n"
      for info in infos:
        if CONTIG in info: 
          important_info += info.strip() + "\n"
          if "_388_" in info:
            print important_info
        elif REAL_ALIGN in info or TOTAL_ALIGN in info:
          which_align = info.split(":")[0].strip()
          lst = info.split(":")[1].split("|")
          gen_align = lst[0].strip().split(" ")
          contig_align = lst[1].strip().split(" ")
          gen_align_begin = int(gen_align[0])
          gen_align_end = int(gen_align[1])
          contig_align_begin = int(contig_align[0])
          contig_align_end = int(contig_align[1])
          first_pos =min( first_pos, gen_align_begin)
            
          if ((gen_align_begin >= begin_pos) and (gen_align_begin <= end_pos) ) or ((gen_align_end >= begin_pos)and(gen_align_end <= end_pos)):
            between_poses = True
          important_info += which_align + ":" + str(gen_align_begin) + " " + str (gen_align_end) + " | " + str(contig_align_begin) + " " + str(contig_align_end) + "\n"
        elif MISASSEMBl in info:
          important_info += info.strip() + "\n"
        #out.write(info)
      if between_poses:
        result_info.append((first_pos, important_info))
        #out.write(important_info)
    infos = [line]
  else:
    infos.append(line)
result_info.sort()
for info in result_info:
  out.write(info[1])
out.write("-----------\n")

out.close()

