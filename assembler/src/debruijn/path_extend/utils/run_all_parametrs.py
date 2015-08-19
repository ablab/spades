import sys
import os

if len(sys.argv) <4:
    print ("<pe_param_file> <out_file> <result_file>")
    exit(1)

out_file = sys.argv[2]
result_file = sys.argv[3]

lst = [0.0, 0.0001, 0.0008, 0.001, 0.002, 0.003, 0.004, 0.005, 0.008, 0.01, 0.015, 0.02,  0.04, 0.05, 0.1, 0.15, 0.2,  0.4, 0.5]
lst1 = [0.0001, 0.0008, 0.001, 0.005, 0.01, 0.05, 0.08, 0.1, 0.2]
lst2 = [ 0.05, 0.08, 0.1,0.15, 0.2, 0.3]
lst1 = [0.06, 0.1, 0.14]
lst1 = [0.001, 0.005]
lst2 = [0.1, 0.13, 0.15, 0.2]
lst1 = [-1.0]
for a in lst:
   # for b in lst:
        pe_params = open(sys.argv[1])
        pe_params_new = open("temp.txt", "w")
        line_index = 1
        for line in pe_params:
            new_line = line.strip() + "\n"
            if line_index == 72:
                new_line = "single_threshold           " + str(a) + "\n"
            #elif line_index == 87:
            #    new_line = "single_threshold           " + str(b) + "\n"
            pe_params_new.write(new_line)
            line_index +=1
        print "change lines"
        pe_params_new.close()
        pe_params.close()
        os.rename("temp.txt", sys.argv[1])
        print "rename file"
        os.system("./run > log_" + str(a) + ".txt")#"_" + str(b) + ".txt")
        print "run program"
        os.system("cp -r " + out_file + " " + result_file + "_" + str(a) + ".fasta")#"_"+str(b) + ".fasta")
        print "copy file"+ out_file + " " + result_file + "_" + str(a) +".fasta"# "_"+str(b) + ".fasta"

