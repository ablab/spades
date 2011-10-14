import sys
import os

GENOME_PATH = '/data/EULER/ablab/input/E.Coli.K12.MG1655/'
GENOME_NAME = 'MG1655-K12'
GENOME_SUFF = '.fasta'
GENOME_FILE = GENOME_PATH + GENOME_NAME + GENOME_SUFF

PAIRED_READ_1 = sys.argv[1]

if (PAIRED_READ_1 == '-h'):
	print("Usage: file1 file2 ouput_dir/ crop_len offset min_IS max_IS format(-f/-q)\n")

PAIRED_READ_2 = sys.argv[2]

query = ' -1 '+ PAIRED_READ_1 +' -2 '+ PAIRED_READ_2
if (PAIRED_READ_1 == '--single'):
	query = ' ' + PAIRED_READ_2

INDEX_DIR = 'index/'

OUTPUT_DIR = sys.argv[3]
OUTPUT_PREF = 'E.Coli.'
OUTPUT_SUFF = '.fastq'

BOWTIE_PATH = 'bowtie/bowtie-0.12.7/'
BT_BUILD = 'bowtie-build'
BT = 'bowtie'

max = 100000
if len(sys.argv) >= 5:
	max = int(sys.argv[4])

offset = 2000000
if len(sys.argv) >= 6:
     	offset = int(sys.argv[5])

minIS = 400
if len(sys.argv) >= 7:
       	minIS = int(sys.argv[6])
maxIS = 600
if len(sys.argv) >= 8:
       	maxIS = int(sys.argv[7])

format = 'q'
if len(sys.argv) >= 9 and sys.argv[8] == '-f':
	format = 'f'

islist = ' -c '
if len(sys.argv) >= 10 and sys.argv[9] == '--list':
        islist = ''



ALL_FILE_SUFF = str(max) +'_' + str(offset) +'_'+ str(minIS) +'_'+str(maxIS)

# crop EColi genome to first 'max' basepairs
output_genome = GENOME_PATH + GENOME_NAME + '.' + ALL_FILE_SUFF + GENOME_SUFF # output genome (cropped)
f1 = open(GENOME_FILE)
f2 = open(output_genome, 'w')

for line in f1:
	f2.write(line)
	break

cnt = -offset
if offset > 0:
	for line in f1:
		if cnt >= 0:
			break
		cnt+=70

cnt = 0
for line in f1:
	if (max - cnt < 70):
		line = line[0:max-cnt]
	f2.write(line);
	cnt += 70
	if cnt >= max:
		break

f1.close()
f2.close()

print('Cropped\n')

# build bowtie index
os.system(BOWTIE_PATH + BT_BUILD +' '+ output_genome + ' ' +  INDEX_DIR + OUTPUT_PREF+'_tmp_'+ALL_FILE_SUFF+ ' >> logBT'+ALL_FILE_SUFF)
print('Index built\n')

# align reads using bowtie
os.system(BOWTIE_PATH+BT +' ' + islist + ' -' + format + ' -I ' +  str(minIS)+' -X '+ str(maxIS) +' '+ INDEX_DIR+OUTPUT_PREF+'_tmp_'+ALL_FILE_SUFF 
+ ' ' + query
+' --al '+ OUTPUT_DIR+OUTPUT_PREF + ALL_FILE_SUFF + OUTPUT_SUFF + ' > '+OUTPUT_DIR+'log'+ALL_FILE_SUFF +' 2> '+ OUTPUT_DIR+'errlog'+ALL_FILE_SUFF)
print('Alligned\n')
# use -X for maxins length


# gzip resulting files
#os.system('gzip ' + OUTPUT_DIR + OUTPUT_PREF + ALL_FILE_SUFF + '_1' + OUTPUT_SUFF)
#os.system('gzip ' + OUTPUT_DIR + OUTPUT_PREF + ALL_FILE_SUFF + '_2' + OUTPUT_SUFF)

print('Gzipped\n')
