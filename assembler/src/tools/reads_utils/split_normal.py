#!/usr/bin/python -O

import sys

inFile = open(sys.argv[1])
file1 = open(sys.argv[2], 'w')
file2 = open(sys.argv[3], 'w')
sfile1 = open(sys.argv[4], 'w')
sfile2 = open(sys.argv[5], 'w')

line = inFile.readline()

if not line:
    sys.exit()

id1 = line.split('/')[0];
str1 = int(line.split('/')[1]);

read1 = '';
line = inFile.readline();
while line and not line.startswith('>'):
    read1 += line.strip()
    line = inFile.readline();


while 1:
    if not line:
        if (str1 == 1):
            sfile1.write(id1 + '/1\n')
            sfile1.write(read1 + '\n')        
        else:
            sfile2.write(id1 + '/2\n')
            sfile2.write(read1 + '\n') 
        break
    
    id2 = line.split('/')[0];
    str2 = int(line.split('/')[1]);

    read2 = '';
    line = inFile.readline();
    while line and not line.startswith('>'):
        read2 += line.strip()
        line = inFile.readline();
    
    if (id1 == id2):
        file1.write(id1 + '/1\n')
        file1.write(read1 + '\n')        
        file2.write(id2 + '/2\n')
        file2.write(read2 + '\n')        

        if not line:
            break

        id1 = line.split('/')[0];
        str1 = int(line.split('/')[1]);

        read1 = '';
        line = inFile.readline();
        while not line.startswith('>'):
            read1 += line.strip()
            line = inFile.readline();

    else:
        if (str1 == 1):
            sfile1.write(id1 + '/1\n')
            sfile1.write(read1 + '\n')        
        else:
            sfile2.write(id1 + '/2\n')
            sfile2.write(read1 + '\n') 
        
        id1 = id2
        str1 = str2
        read1 = read2       

inFile.close()
file1.close()
file2.close()
sfile1.close()
sfile2.close()

