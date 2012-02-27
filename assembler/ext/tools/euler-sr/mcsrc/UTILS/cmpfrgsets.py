#! /usr/bin/env python


import sys
import string
import re
import Numeric

def Align(s1, s2):
    s = Numeric.zeros((len(s1)+1, len(s2)+1))
    p = Numeric.zeros((len(s1)+1, len(s2)+1))
    for i in range(1,len(s1)+1):
        for j in range(1,len(s2)+1):
            if (s1[i-1] == s2[j-1]):
                s[i,j] = s[i-1,j-1]
                p[i,j] = 0
            else:
                if(s[i-1,j] < s[i,j-1]):
                    s[i,j] = s[i-1,j] + 1
                    p[i,j] = 1
                else:
                    s[i,j] = s[i,j-1]  + 1
                    p[i,j] = 2
    return s[len(s1),len(s2)]

            
def KAlign(s1,s2,k):
    while (len(s1) < len(s2)):
        s1 = s1 + "-"
    while (len(s2) < len(s1)):
        s2 = s2 + "|"
    middle = k/2
    s = Numeric.zeros((len(s2)+1,k))
    end = len(s1)
    p = 0
    for p in range(1,len(s2)):
        for q in range (0,k):
            if (p + q - middle < 0):
                continue
            elif (p + q - middle >= end):
                continue
            else:
                i = p + q - middle
                j = p
                if (s1[i] == s2[j]):
                    s[p,q] = s[p-1,q]
                else:
                    # not equal.  Handle in order of boundary
                    #conditions.
                    # left boundary condition
                    if (q == 0):
                        s[p,q] = s[p-1,q+1] + 1 # score from above
                    # middle boundary condition
                    elif (q < k-1):
                        if (s[p,q-1] < s[p-1,q+1]):
                            s[p,q] = s[p,q-1] + 1  # come from left
                        else:
                            s[p,q] = s[p-1,q+1] + 1 # come from above
                    else:
                        # right boundary condition
                        s[p,q] = s[p,q-1] + 1
    return s[p,k/2]/2
            
            
            


if (len(sys.argv) < 2): 
    print "usage: cmpfrgsets  file1.seq file2.seq";
    sys.exit(2)
else:
    file1Name = sys.argv[1]
    file2Name = sys.argv[2]

f1 = open(file1Name);
f2 = open(file2Name);

l1 = ""
l2 = ""
line = 1;

net = 0
tl1 = ""
tl2 = ""
frag = 0
while ((tl1 != None) and (tl2 != None)):
    tl1 = f1.readline().strip()
    tl2 = f2.readline().strip()
    while ((tl1 != None) and (tl1[0] != '>')):
        l1 = l1 + tl1
        tl1 = f1.readline().strip()
        if (tl1 == ""):
            tl1 = None
    while ((tl2 != None) and (tl2[0] != '>')):
        l2 = l2 + tl2
        tl2 = f2.readline().strip()
        if (tl2 == ""):
            tl2 = None

    net = net + KAlign(l1,l2,9)

    frag = frag + 1
        
    l1 = ""
    l2 = ""

print "Net aligned differences: ", net
