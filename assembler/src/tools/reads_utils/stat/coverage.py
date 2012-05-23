#!/usr/bin/python -O

#Calculate coverage from raw file

import sys

def coverage(in_filename, out_filename, maxLen, bar, k):
    
    inFile = open(in_filename)
    outFile = open(out_filename, 'w')
    
    hist = [0 for i in range(maxLen	+ 1)]

    for line in inFile:
	    stpos = int(line.split()[0])
	    for i in range(0, int(line.split()[1]) - k + 1):
		    cpos = stpos + i
		    if cpos <= maxLen:
			    hist[cpos] += 1

    covered = 0.0
    for i in range(0,maxLen + 1):
	    if (hist[i] > 0):
		    covered += 1.0

    # print("Coverage: " + str(covered/maxLen) + "\n")

    newHist = [0 for i in range((maxLen + 1) / bar + 2)]

    for i in range(maxLen + 1):
	    newHist[int(i/bar)] += hist[i]

    for i in range((maxLen + 1) / bar + 1):
	    outFile.write(str(i) + ' ' + str(newHist[i] / bar) + '\n')

    inFile.close()
    outFile.close()

    return covered/(maxLen + 1)


def analyze_gaps(in_filename, out_filename, reference, out_ref, k)
    inFile = open(in_filename)
    outFile = open(out_filename, 'w')

    gaps100 = 0
    gaps500 = 0
    gaps1000 = 0

    current = 0

    #[start,end)
    gaps = []
    chunks = []

    line = inFile.readline()
    while line:
	    cov = int(line.split()[1])

            end = current
            while cov == 0 and line:
                   end += 1 
                   line = inFile.readline()
                   if line:
                          cov = int(line.split()[1])

            if end != current:
                   gaps.append((current, end));

            current = end
            while cov > 0 and line:
                   end += 1 
                   line = inFile.readline()
                   if line:
                          cov = int(line.split()[1])

            if end != current:
                   chunks.append((current, end + k - 1));

            current = end

def main():

    if len(sys.argv) < 5:
	    print("Usage: <coverage file> <output> <genome length> <bar width> [k = 1]");
	    exit(0)

    k = 1
    if len(sys.argv) > 5:
	k = int(sys.argv[5])

    cov = coverage(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), k)

    print("Coverage: " + str(cov) + "\n")

if __name__ == '__main__':
    main()
