#!/usr/bin/python -O

#Calculate coverage from raw file

import sys
import os
import datetime

profile = []
insertions = {}
def read_genome(filename):
    res_seq = []

    seq =''
    for line in open(filename):
            if line[0] == '>':
    	    	    res_seq.append(line[1:])
            else:
    	    	    seq += line.strip()
    res_seq.append(seq)
    return res_seq


def write_fasta(data, filename):
    outFile = open(filename, 'w')
    #print data[0][0];
    for seq in data:
#            print(seq[0])
            outFile.write('>' + seq[0] );
            i = 0
            while i < len(seq[1]):
                    outFile.write(seq[1][i:i+60] + '\n')
                    i += 60

    outFile.close()

def vote_insertions(position):
    global insertions;
    arr = insertions[position]
    lengths = {}
    for strn in arr:
        if len(strn) not in lengths:
            lengths[len(strn)] = 0;
        lengths[len(strn)] += 1;
    opt_len = 0;
    max = 0;
    for l in lengths:
        if lengths[l] > max:
            max = lengths[l]
            opt_len = l;
    print "insertion; length =" + str(opt_len)
    ins_profile = []
    for i in range (0, opt_len):
        ins_profile.append({})
        for j in ('A','C','T','G','N'):
            ins_profile[i][j] = 0;
    for strn in arr:
#Others are ignored
        if len(strn) == opt_len:
            for i in range (0, opt_len):
                ins_profile[i][strn[i]]+=1;
    best_ins = ''
    for i in range (0, opt_len):
        next_nucl = 'N';
        for j in ('A','C','T','G'):
            if ins_profile[i][j] > ins_profile [i][next_nucl]:
                next_nucl = j;
        best_ins += next_nucl;
    return best_ins


def process_read(cigar, aligned, position, l):
    global profile
    global insertions
    l_read = len(aligned);

#    if position > 36470:
#        print aligned+ " "+ cigar +" " + str(position)
#    if 'I' in cigar or 'D' in cigar or '*' in cigar:
    if '*' in cigar:
        return 0;

    shift = 0;
    nums = []
    operations = []
    tms = "";
    aligned_length = 0;
    for i in range(0, len(cigar)):

        if cigar[i].isdigit():
             tms += cigar[i]
        else:
            nums.append(int(tms));
            operations.append(cigar[i]);
            if cigar[i]=='M':
                aligned_length += int(tms)
            
	    if (cigar[i] =='D' or cigar[i] == 'I') and int(tms) > 5:
#TODO: dirty hack, we do not fix big indels now
		return 0;
	    tms = ''
#we do not need short reads aligned near gaps
    if aligned_length < 40 and position > 50 and l - position > 50:
        return 0;
    state_pos = 0;
    shift = 0;
    skipped = 0;
#    if position < 150:
#        print aligned+ " "+ cigar +" " + str(position)
    insertion_string = ''
    for i in range(0, l_read):
    #            print aligned[i];
    #            print profile[i+position - 1]
#        if position == 698:
#            print "inserting" + str(i) +" " + str(position) + " " + str(skipped)+ " " + str(l_read) + " " + cigar +" " + str(l) + " " + insertion_string + "  " + operations[state_pos]

        if i + position - skipped >= l:
#            print str(i) +" " + str(position) + " " + str(skipped)+ " " + str(l_read) + " " + cigar +" " + str(l) + " " + insertion_string + "  " + operations[state_pos]
#            print "read going out of contig"
            break;

        if shift + nums[state_pos] <= i:
            shift += nums[state_pos]
            state_pos += 1
#        if position == 698:
#            print operations[state_pos]
        if insertion_string != '' and operations[state_pos] != 'I':
#            print "inserting" + str(i) +" " + str(position) + " " + str(skipped)+ " " + str(l_read) + " " + cigar +" " + str(l) + " " + insertion_string + "  " + operations[state_pos]
            ind = i+position-skipped - 1;
            if ind not in insertions:
                insertions[ind]= []
            insertions[ind].append(insertion_string)
            insertion_string = ''
        if operations[state_pos] == 'M':
            if i +  position - skipped < l:
                profile[i+position -skipped][aligned[i]] += 1
        else:
            if  operations[state_pos] in {'S', 'I', 'H'}:
                if operations[state_pos] == 'I':
#                    print "ddd" + str(i) +" " + str(position) + " " + str(skipped)+ " " + str(l_read) + " " + cigar +" " + str(l) + " " + insertion_string + "  " + operations[state_pos]
                    if insertion_string == '':
                        profile[i+position -skipped - 1 ]['I'] += 1
                    insertion_string += (aligned[i])
                skipped += 1;
            elif operations[state_pos] == 'D':
#               print str(i) +" " + str(position) + " " + str(skipped)+ " " + str(l_read) + " " + cigar +" " + str(l) + " " + insertion_string + "  " + operations[state_pos]
                profile[i+position - skipped]['D'] += 1
                skipped -=1;
                shift -= 1;
    if insertion_string != '' and operations[state_pos] != 'I':
#        print "inserting" + str(i) +" " + str(position) + " " + str(skipped)+ " " + str(l_read) + " " + cigar +" " + str(l) + " " + insertion_string + "  " + operations[state_pos]

        ind = i+position-skipped - 1;
        if ind not in insertions:
            insertions[ind]= []
        insertions[ind].append(insertion_string)
    return 1
def split_fasta(filename):
#    for
    return 32

def main():
    global profile
    global insertions
    if len(sys.argv) < 2:
	    print("Usage: dir with .pair.sam and fasta files")
	    exit(0)

    replaced = 0;
    inserted = 0;
    deleted = 0;
    now = datetime.datetime.now()
    res_directory = "corrector.output." + now.strftime("%Y.%m.%d_%H.%M.%S")+"/";
    if not os.path.exists(res_directory):
        os.makedirs(res_directory)
#    refinedFileName = res_directory + sys.argv[2].split('/')[-1].split('.')[0] + '.ref.fasta';
    dir = sys.argv[1];

    filelist = [os.path.abspath(os.path.join(dir, i)) for i in os.listdir(dir) if os.path.isfile(os.path.join(dir, i))]
    for contig_file in filelist:
        f_name = contig_file.split('/')[-1];
        f_arr = f_name.split('.');
 #       print contig_file + "  is contig_file";
 #       print os.path.join("/".join(contig_file.split('/')[:-1]), f_arr[0]+".pair.sam")

        if len(f_arr) == 2 and f_arr[1] == "fa" and len(f_arr[0]) < 16 and os.path.exists(os.path.join("/".join(contig_file.split('/')[:-1]), f_arr[0]+".pair.sam")):

            samfilename = os.path.join("/".join(contig_file.split('/')[:-1]), f_arr[0]+".pair.sam");

            samFile = open(samfilename, 'r');
            fasta_contig = read_genome(contig_file);
            print "processing " + str(contig_file) + ", contig length:" + str(len(fasta_contig[1]));
            contig = fasta_contig[1].upper()
#            profile = []
            total_reads = 0;
            indelled_reads = 0;
#            print samfilename
            refinedFileName = res_directory + samfilename.split('/')[-1].split('.')[0] + '.ref.fasta';
    #    samFile = open( sys.argv[1], 'r');
# accurate!
            cont_num = contig_file.split('/')[-1].split('.')[0];
#    fasta_contig = read_genome(sys.argv[2]);



            l = len(contig)
#            print l;
#            print contig
            rescontig = ""
            for i in range (0, l):
                profile.append( {} );
                for j in ('A','C','G','T','N','I','D'):
                    profile[i][j] = 0;
            insertions = {}
#TODO: estimation on insert size, need to be replaced            
            insert_size_est = 500
            for line in samFile:
                arr = line.split();
                if arr[2].split('_')[1] != cont_num:
                    continue
        #        print line;
                cigar = arr[5]
                aligned = arr[9];
                position = int(arr[3]) - 1
		mate = arr[6];
		if mate != '=' and mate != '*' and position > insert_size_est and position < l - insert_size_est - 100:
		    continue;
        #        print position;
        #        print l;
                indelled_reads += 1 - process_read(cigar, aligned, position, l)
                total_reads += 1;
            if len(insertions) < 50:
                print insertions
            else:
                print "insertions very big, most popular:"
                for element in insertions:
                    if len(insertions[element]) > 10:
                        print str(element) + str(insertions[element])
                        print "profile here: " + str(profile[element])
            for i in range (0, l):
        #        print profile[i];
                tj = contig[i]
                tmp =''
                for count in range(0,2):
                    for j in ('A','C','G','T','N', 'I', 'D'):
                        if profile[i][tj] < profile[i][j] or (j == 'I' and profile[i][tj] < 1.5 * profile[i][j] and profile[i][j] != 1):
                            tj = j
            #                rescontig[i] = j

                    if tj != contig[i] :
                        if tj =='I' or tj == 'D' :
                            print "there was in-del"
                        else:
                            print profile[i]
                            print "changing " + contig[i] + " to " + tj + " on position " + str(i+1) + '\n'
                            replaced += 1;
                    if tj in ('A','C','T','G','N'):
                        rescontig += tj
                        break;
                    elif tj == 'D':
                        print "skipping deletion"+ " on position " + str(i+1) + '\n'
                        deleted += 1;
                        break;
                    elif tj == 'I':
                        tmp = vote_insertions(i);

                        profile[i]['I'] = 0;
                if tmp != '':
                    rescontig += tmp;
                    inserted += len (tmp)

                #            print(fasta_contig[0]);
            res_fasta = [[]]
            #res_fasta.append([])
            res_fasta[0].append(fasta_contig[0]);
            res_fasta[0].append(rescontig)
            write_fasta(res_fasta, refinedFileName)
        #    refinedFile.write(rescontig);
 #           nonFasta.write(rescontig);
            profile = []
            print "total/indelled reads:" + str(total_reads) + '/' + str(indelled_reads)
    print "TOTAL replaced: "+ str(replaced) + " inserted: "+ str(inserted) + " deleted: " + str(deleted);

if __name__ == '__main__':
    main()
