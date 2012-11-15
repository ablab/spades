#!/usr/bin/python -O

#Calculate coverage from raw file

import sys
import os
import datetime
import io
import cProfile
import getopt
#from joblib import Parallel, delayed

profile = []
insertions = {}
def read_genome(filename):
    res_seq = []

    seq =''
    for line in io.open(filename):
            if line[0] == '>':
    	    	    res_seq.append(line[1:])
            else:
    	    	    seq += line.strip()
    res_seq.append(seq)
    return res_seq

def read_contigs(filename):
    res_seq = {}

    seq =''
    cont_name = ''
    for line in io.open(filename):
        if line[0] == '>':
            if cont_name != '':
                res_seq[cont_name] = seq
            seq = ''
            cont_name = line[1:].strip()
        else:
            seq += line.strip()
    if cont_name != '':
        res_seq[cont_name] = seq
    return res_seq


def write_fasta(data, filename):
    outFile = io.open(filename, 'w')
    for seq in data:
            outFile.write('>' + seq[0] );
            i = 0
            while i < len(seq[1]):
                    outFile.write(seq[1][i:i+60] + '\n')
#                    print seq[1][i:i+60]
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


def process_read(cigar, aligned, position, l, mate):
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
    if aligned_length < 40 and position > 50 and l - position > 50 and mate == 1:
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
                profile[i+position -skipped][aligned[i]] += mate
        else:
            if  operations[state_pos] in {'S', 'I', 'H'}:
                if operations[state_pos] == 'I':
#                    print "ddd" + str(i) +" " + str(position) + " " + str(skipped)+ " " + str(l_read) + " " + cigar +" " + str(l) + " " + insertion_string + "  " + operations[state_pos]
                    if insertion_string == '':
                        profile[i+position -skipped - 1 ]['I'] += mate
                    insertion_string += (aligned[i])
                skipped += 1;
            elif operations[state_pos] == 'D':
#               print str(i) +" " + str(position) + " " + str(skipped)+ " " + str(l_read) + " " + cigar +" " + str(l) + " " + insertion_string + "  " + operations[state_pos]
                profile[i+position - skipped]['D'] += mate
                skipped -=1;
                shift -= 1;
    if insertion_string != '' and operations[state_pos] != 'I':
#        print "inserting" + str(i) +" " + str(position) + " " + str(skipped)+ " " + str(l_read) + " " + cigar +" " + str(l) + " " + insertion_string + "  " + operations[state_pos]

        ind = i+position-skipped - 1;
        if ind not in insertions:
            insertions[ind]= []
        insertions[ind].append(insertion_string)
    return 1
def split_sam(filename, tmpdir):

    inFile = io.open(filename, buffering=131072)
    separate_sams ={}
    print("file io.opened")
    for line in inFile:
        arr = line.split('\t')
        if len(arr) > 5:
            if arr[2] in separate_sams:
                separate_sams[arr[2]].write(line)
        else:
            contig = arr[1].split(':')
            if contig[0] == 'SN':
                samfilename = tmpdir + '/' +contig[1].split('_')[1] + '.pair.sam'
#                print contig[1] + " " + samfilename;
                separate_sams[contig[1]] = io.open(samfilename, 'w', buffering=131072)

    for file_name in separate_sams:
        separate_sams[file_name].close()

    return 0
def split_contigs(filename, tmpdir):
    ref_seq = read_contigs(filename);
    for contig_desc in ref_seq:
        tfilename = tmpdir + '/' + contig_desc.split('_')[1] +'.fasta'

        write_fasta([[contig_desc, ref_seq[contig_desc]]], tfilename)

    return 0

def run_bwa():
# align with BWA
    global config;
    os.system("bwa index -a is " + config.contigs + " 2")
#    (contigs_name, path, suf) = fileparse(config.contigs)

    os.system("bwa aln  "+ contigs_name +" " +  config.reads1 + " -t 4 -O 7 -E 2 -k 3 -n 0.08 -q 15")
    os.system("bwa aln  "+ contigs_name +" " +  config.reads2 + " -t 4 -O 7 -E 2 -k 3 -n 0.08 -q 15")
    os.system("bwa sampe "+ config.contigs + " " + config.sai1 + " " + config.reads1 +" " + config.sai2 + " " + config.reads2 + "> tmp.sam 2 > isize.txt")
#    my ($sai1, $sai2, $sam) = ("reads1_aln.sai", "reads2_aln.sai", "reads_aln.sam");
#    my $cmd;
#    $cmd = "$bwa index -a is $contigs 2>/dev/null";
#    system($cmd);
#
#    my ($max_ge, $max_go, $k, $n) = (100, 3, 3, 0.08);
#    $bwa_aln1_cmd  = "$bwa aln $contigs_name $reads1 -t $num_threads -O 7 -E 2 -k $k -n $n -q 15";
#    $bwa_aln2_cmd  = "$bwa aln $contigs_name $reads2 -t $num_threads -O 7 -E 2 -k $k -n $n -q 15";
#}
#    $cmd = $bwa_aln1_cmd . " > $sai1";
#    system($cmd);
#    $cmd = $bwa_aln2_cmd . " > $sai2";
#    system($cmd);
#
#    print "\nGenerating paired-end alignments & estimating insert-size...\n";
#    my $is = "isize.txt";
#    my $bwa_sampe_cmd = "$bwa sampe $contigs $sai1 $sai2 $reads1 $reads2 > $sam 2>$is";
    return 0;
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
    tmpdir = res_directory +'tmp';
#    tmpdir = 'tmp'
#    os.makedirs(tmpdir)
#    split_contigs(sys.argv[2], tmpdir)
#    print("contigs splitted, starting splitting samfile");
#    split_sam(sys.argv[1], tmpdir)
#    return 0
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

            samFile = io.open(samfilename, 'r');
            fasta_contig = read_genome(contig_file);
            print "processing " + str(contig_file) + ", contig length:" + str(len(fasta_contig[1]));
            contig = fasta_contig[1].upper()
#            profile = []
            total_reads = 0;
            indelled_reads = 0;
#            print samfilename
            refinedFileName = res_directory + samfilename.split('/')[-1].split('.')[0] + '.ref.fasta';
    #    samFile = io.open( sys.argv[1], 'r');
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

                mate_el = arr[6];
                if abs(position - 188264) < 100:
                    print arr[0] + " "+ arr[3] +" "+ mate_el
                if mate_el != '=' and mate_el != '*' and position > insert_size_est and position < l - insert_size_est - 100:
                    if abs(position - 188264) < 100:
                        print "skipping"
                    continue;
        #        print position;
        #        print l;

                if mate_el == '=' and (int(arr[1]) & 8) == 0:
                    mate = 2;
                else:
                    mate = 1;
                indelled_reads += 1 - process_read(cigar, aligned, position, l, mate)
                total_reads += 1;
            if len(insertions) < 50:
                print insertions
            else:
                print "insertions very big, most popular:"
                for element in insertions:
                    if len(insertions[element]) > 30:
                        print str(element) + str(insertions[element])
                        print "profile here: " + str(profile[element])
            for i in range (0, l):
                if i in {188260, 188264, 188268, 188272, 188276, 188248}:
                    print "FFFFFFUUUU "+ str(i + 1)
                    print(profile[i])
        #        print profile[i];
                tj = contig[i]
                tmp =''
                for count in range(0,2):
                    for j in ('A','C','G','T','N', 'I', 'D'):
                        if profile[i][tj] < profile[i][j] or (j == 'I' and profile[i][tj] < 1.5 * profile[i][j] and profile[i][j] > 2):
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
