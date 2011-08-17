import sys
import os

DATA_PATH = '../../data/input'
BOWTIE_PATH = 'bowtie/bowtie-0.12.7'

max = 100000
if len(sys.argv) >= 2:
	max = int(sys.argv[1])

# crop EColi genome to first 'max' basepairs
filename1 = DATA_PATH + '/MG1655-K12.fasta' # input genome
filename2 = 'MG1655-K12.' + str(max) + '.fasta' # output genome (cropped)
f1 = open(filename1)
f2 = open(filename2, 'w')
cnt = -70; # because of the '>' in first line
for line in f1:
	if (max - cnt < 70):
		line = line[0:max-cnt]
	f2.write(line);
	cnt += 70
	if cnt >= max:
		break
f1.close()
f2.close()

# build bowtie index
os.system(BOWTIE_PATH + '/bowtie-build ' + filename2 + ' bowtie/ecoli.'+str(max))

# align reads using bowtie
os.system(BOWTIE_PATH + '/bowtie bowtie/ecoli.'+str(max) + ' -1 ' + DATA_PATH + '/EAS20_8_Quake/s_6_1.cor.fastq -2 ' + DATA_PATH + '/EAS20_8_Quake/s_6_2.cor.fastq --al s_6.first' + str(max) + '.fastq')
# use -X for maxins length

# gzip resulting files
os.system('gzip s_6.first' + str(max) + '_1.fastq')
os.system('gzip s_6.first' + str(max) + '_2.fastq')
