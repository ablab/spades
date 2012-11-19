import sys

if len(sys.argv) != 4:
  print "<first quast loose genes file> <second quast loose genes file> <output>"
  exit(1)

first = open(sys.argv[1], 'r')
second = open(sys.argv[2], 'r')
out = open(sys.argv[3], 'w')

PART = 'part'
TOTAL = 'total'
spliter = '\t'
def is_part(line):
  if PART in line:
    return True
  return False

first_part = dict()
second_part = dict()
first_total = dict()
second_total = dict()

def find_types(quast_file, part, total):
  for line in quast_file:
    info = line.split(spliter)
    loose_id = int(info[0])
    if is_part(line):
      part[loose_id] = line.strip()
    else:
      total[loose_id] = line.strip()

find_types(first, first_part, first_total)
find_types(second, second_part, second_total)

out.write("First loose\n")
for loose_id, line in first_part.items():
  if loose_id not in second_part and loose_id not in second_total:
    out.write(line + " " + PART + "\n")
for loose_id, line in first_total.items():
  if loose_id not in second_total:
    if loose_id not in second_part: 
      out.write(line + " " + TOTAL + "\n")
    else:
      out.write(line + " "  + TOTAL + " " + PART + "\n")
out.close()


