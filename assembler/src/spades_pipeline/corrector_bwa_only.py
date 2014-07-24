#!/usr/bin/python -O

############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#Calculate coverage from raw file
import glob
import gzip
import sys
import os
import datetime
import getopt
from site import addsitedir
import logging
import shutil

from math import pow
from support import universal_sys_call, error
import options_storage

#profile = []
#insertions = {}
config = {}
#total_contigs


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


def read_contigs(filename):
    res_seq = {}

    seq =''
    cont_name = ''
    for line in open(filename):
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
    outFile = open(filename, 'w')
    for seq in data:
            outFile.write('>' + seq[0].strip() + '\n' )
            i = 0
            while i < len(seq[1]):
                    outFile.write(seq[1][i:i+60] + '\n')
                    i += 60
    outFile.close()


def process_read(arr, l,m, profile, insertions):
#    global profile
#    global insertions
    cigar = arr[5]
    aligned = arr[9]
    qual_aligned = arr[10]
    position = int(arr[3]) - 1
    map_q = int(arr[4])
    mate = m
    l_read = len(aligned)

    if '*' in cigar:
        return 0

    if config["use_quality"]:
        mate *=  (1 - pow(10, map_q / -10))
    nums = []
    operations = []
    tms = ""
    aligned_length = 0
    for i in range(0, len(cigar)):

        if cigar[i].isdigit():
             tms += cigar[i]
        else:
            nums.append(int(tms))
            operations.append(cigar[i])
            if cigar[i]=='M':
                aligned_length += int(tms)
            tms = ''
#we do not need short reads aligned near gaps
    if aligned_length < min(l_read* 0.4, 40) and position > l_read / 2 and l - position > l_read / 2 and mate == 1:
        return 0
    state_pos = 0
    shift = 0
    skipped = 0
    deleted = 0
#    if position < 150:
#        print aligned+ " "+ cigar +" " + str(position)
    insertion_string = ''
    int_A = ord('A') - 1
    for i in range(0, l_read):
    #            print aligned[i]
    #            print profile[i+position - 1]
        if i + position - skipped >= l:
#            print str(i) +" " + str(position) + " " + str(skipped)+ " " + str(l_read) + " " + cigar +" " + str(l) + " " + insertion_string + "  " + operations[state_pos]
#            print "read going out of contig"
            break

        if shift + nums[state_pos] <= i:
            shift += nums[state_pos]
            state_pos += 1
#        if position == 698:
#            print operations[state_pos]
        if insertion_string != '' and operations[state_pos] != 'I':
#            print "inserting" + str(i) +" " + str(position) + " " + str(skipped)+ " " + str(l_read) + " " + cigar +" " + str(l) + " " + insertion_string + "  " + operations[state_pos]
            ind = i + position - skipped - 1
            if ind not in insertions:
                insertions[ind]= []
            insertions[ind].append(insertion_string)
            insertion_string = ''
        if operations[state_pos] == 'M':
            if i +  position - skipped < l:
                t_mate = mate
                if config["use_quality"]:
                    k = pow(10, (ord(qual_aligned[i - deleted]) - int_A ) / -10)
                    if qual_aligned[i - deleted] == 'B':
                        t_mate *= 0.2
                    else:
                        t_mate *= (1 - k / (k + 1))
                profile[i + position - skipped][aligned[i - deleted]] += t_mate
        else:
            if  operations[state_pos] in ('S', 'I', 'H'):
                if operations[state_pos] == 'I':
                    if insertion_string == '':
                        profile[i + position - skipped - 1]['I'] += mate
                    insertion_string += (aligned[i])
                skipped += 1
            elif operations[state_pos] == 'D':
#               print str(i) +" " + str(position) + " " + str(skipped)+ " " + str(l_read) + " " + cigar +" " + str(l) + " " + insertion_string + "  " + operations[state_pos]
#                print str(i) + " " + str(i+position - skipped) + "    deletion!"
                profile[i + position - skipped]['D'] += mate
#                skipped =1
#                shift -= 1
                deleted += 1
    if insertion_string != '' and operations[state_pos] != 'I':
#        print "inserting" + str(i) +" " + str(position) + " " + str(skipped)+ " " + str(l_read) + " " + cigar +" " + str(l) + " " + insertion_string + "  " + operations[state_pos]

        ind = i + position - skipped - 1
        if ind not in insertions:
            insertions[ind]= []
        insertions[ind].append(insertion_string)
    return 1


def drop_cash(cashed_sams, separate_sams, log, last):
    log.info("dropping cashed lines...")
    limit = 0
    if not last:
        limit = 100
    for contig_info in cashed_sams:
        file_name = separate_sams[contig_info]
        if len(cashed_sams[contig_info]) > limit:
            file = open(file_name, "a")
            for line in cashed_sams[contig_info]:
                file.write(line)
            file.close()
            cashed_sams[contig_info] = []


def split_sam(filename, tmpdir, log):
    #global total_contigs
    separate_sams ={}
    mult_aligned = {}
    cashed_sams = {}
    log.info("file " + filename + " opened")
    read_num = 0

#    import subprocess
#    ulimit = subprocess.check_output("ulimit -n", shell=True)
#    print "ulimit :", ulimit
#    ulimit = int(ulimit)
    need_to_cashe = True
#TODO: whether multiple aligned helps?
#    needed_ulimit = total_contigs * 2 + 20
#    if (ulimit < needed_ulimit):
#        print ("Limit for simultaniously opened file is too small. limit: " + str(ulimit) + "; needed: " + str(needed_ulimit))
#        print ("Increasing it will speedup this stage for a while")
#        print ("For most systems (if you are a sudoer) you can fix it by adding \"*  soft    nofile    "+ str(needed_ulimit + 20) + "\" to /etc/security/limits.conf")
#        need_to_cashe = True

    paired_read = []
    processed_reads = 0
    max_cashed_lines = 1000000
    for line in open(filename, 'r'):
        arr = line.split()
        processed_reads += 1
        if processed_reads % max_cashed_lines == 0 and need_to_cashe:
            log.info(str(processed_reads))
            drop_cash(cashed_sams , separate_sams, log, False)
        if len(arr) > 5:
            read_num += 1
            paired_read.append(line.strip())
            if read_num == 2:
                unique = {}
                all_contigs = {}
                for j in range(0, 2):
                    tags = paired_read[j].split()[11:]
                    parsed_tags = {}
                    unique_fl = paired_read[j].split()[2]
                    all_contigs[unique_fl] = 1
                    AS = 0
                    XS = -1000
                    if paired_read[j].split()[5] == "*":
                        unique_fl = ''
                    else:
                        for tag in tags:
                            tr = tag.split(':')
                            tag_name = tr[0]
                            tag_val = tr[2]
                            if (tag_name == "X0" and int(tag_val) > 1):
                                unique_fl = ''
                            elif (tag_name == "AS"):
                                AS = tag_val
                            elif (tag_name == "XS"):
                                XS = tag_val
                    if unique_fl != '' and AS > XS:
                        unique[unique_fl] = 1

                for cont_name in all_contigs :
                    if cont_name in separate_sams:
                        for j in range(0,2):
                            if (not need_to_cashe):
                                separate_sams[cont_name].write(paired_read[j]+ '\n')
#                                if cont_name in unique:
#                                    separate_sams[cont_name].write(paired_read[j]+ '\n')
#                                else:
#                                    mult_aligned[cont_name].write(paired_read[j]+ '\n')
                            else:
                                cashed_sams[cont_name].append(paired_read[j] + '\n')

                read_num = 0
                paired_read = []
#            if arr[2] in separate_sams:
#                separate_sams[arr[2]].write(line)

        else:
            contig = arr[1].split(':')
            if contig[0] == 'SN':
                samfilename = os.path.join(tmpdir, contig[1].split('.')[0] + '.pair.sam')
                multalignedfilename = os.path.join(tmpdir, contig[1].split('.')[0] + '.multiple.sam')
                #samfilename = tmpdir + '/' +contig[1].split('_')[1] + '.pair.sam'
                #multalignedfilename = tmpdir + '/' +contig[1].split('_')[1] + '.multiple.sam'
#               print contig[1] + " " + samfilename
                if need_to_cashe:
                    separate_sams[contig[1]] = samfilename
                    cashed_sams[contig[1]] = []
                    cashed_sams[contig[1]].append(line)

                else:
                    separate_sams[contig[1]] = open(samfilename, 'w')
                    mult_aligned[contig[1]] = open(multalignedfilename , 'w')
                    separate_sams[contig[1]].write(line)

    if need_to_cashe:
        drop_cash(cashed_sams , separate_sams, log, True)
        cashed_sams = None
    else:
        for file_name in separate_sams:
            separate_sams[file_name].close()

    return 0


def split_reads(filename, left, right):
    #left = os.path.join(config["work_dir"], os.path.splitext(os.path.basename(filename))[0] + "_1.fastq")
    #right = os.path.join(config["work_dir"], os.path.splitext(os.path.basename(filename))[0] + "_2.fastq")
    leftFile = open(left, "a")
    rightFile = open(right, "a")
    count = 0
    state = True
    for line in open(filename, 'r'):
        if state:
            leftFile.write(line)
        else:
            rightFile.write(line)
        count += 1
        if count == 4:
            count = 0
            state = not state
    leftFile.close()
    rightFile.close()
    #config["reads1"] = left
    #config["reads2"] = right


def split_contigs(filename, tmpdir):
    #global total_contigs
    ref_seq = read_contigs(filename)
    for contig_desc in ref_seq:
        tfilename = os.path.join(tmpdir, contig_desc.split('.')[0] + '.fasta')
        #tfilename = tmpdir + '/' + contig_desc.split('_')[1] +'.fasta'
        #total_contigs += 1
        write_fasta([[contig_desc, ref_seq[contig_desc]]], tfilename)


def usage():
    sys.stderr.write('Mismatch Corrector. Simple post processing tool\n')
    sys.stderr.write('Usage: python' + str(sys.argv[0]) + '[options] -1 left_reads -2 right_reads -c contigs\n')
    sys.stderr.write('Or: python' + str(sys.argv[0]) + '[options] --12 mixed_reads -c contigs\n')
    sys.stderr.write('Or: python' + str(sys.argv[0]) + '[options] -s sam_file -c contigs\n')
    sys.stderr.write('Options:\n')
    sys.stderr.write('-t/--threads  <int>   threads number\n')
    sys.stderr.write('-o/--output-dir   <dir_name>  directory to store results\n')
    sys.stderr.write('-m/--mate-weight  <int>   weight for paired reads aligned properly. By default, equal to single reads weight (=1)\n')
    sys.stderr.write('--bwa <path>   path to bwa tool. Required if bwa is not in PATH\n')
    sys.stderr.write('--bowtie2    <path>  path to bowtie2 tool. Can be used instead of bwa. It is faster but provides a bit worse results\n')
    sys.stderr.write('  --use-quality use quality values as probabilities \n')
    sys.stderr.write('  --debug   save all intermediate files \n')
    sys.stderr.write('  --use-multiple-aligned  use paired reads with multiple alignment\n')
    sys.stderr.write('  --skip-masked   do not correct single \'N\' in contigs unless significant read support provided\n')
    sys.stderr.write('  --insert-size <int> estimation on insert size\n')
    sys.stderr.flush()


def run_aligner(log):
    global config

    #    (contigs_name, path, suf) = fileparse(config.contigs)
    if not "contigs" in config:
        error("contigs were not specified!", log)
    if (not "reads1" in config or not "reads2" in config) and (not "reads_mixed" in config):
        error("reads were not specified!", log)

    # splitting and merging all provided reads files to get only two files with left and right reads respectively
    left = os.path.join(config["work_dir"], "reads_left.fastq")
    right = os.path.join(config["work_dir"], "reads_right.fastq")
    if "reads_mixed" in config:
        for mixed_reads_file in config["reads_mixed"]:
            split_reads(mixed_reads_file, left, right)
    if "reads1" in config and (len(config["reads1"]) > 1 or os.path.isfile(left)):
        for reads_file in config["reads1"]:
            if os.path.splitext(reads_file)[1] == '.gz' or (reads_file in options_storage.dict_of_prefixes and
                                                           options_storage.dict_of_prefixes[reads_file].endswith('.gz')):
                fdscr = gzip.open(reads_file, 'r')
            else:
                fdscr = open(reads_file, 'r')
            shutil.copyfileobj(fdscr, open(left, 'a'))
    if "reads2" in config and (len(config["reads2"]) > 1 or os.path.isfile(right)):
        for reads_file in config["reads2"]:
            if os.path.splitext(reads_file)[1] == '.gz' or (reads_file in options_storage.dict_of_prefixes and
                                                           options_storage.dict_of_prefixes[reads_file].endswith('.gz')):
                fdscr = gzip.open(reads_file, 'r')
            else:
                fdscr = open(reads_file, 'r')
            shutil.copyfileobj(fdscr, open(right, 'a'))

    if os.path.isfile(left) and os.path.isfile(right):
        config["reads1"] = left
        config["reads2"] = right
    else:
        config["reads1"] = config["reads1"][0]
        config["reads2"] = config["reads2"][0]

    shutil.copy(config["contigs"], config['work_dir'])
    #os.system("cp " + config["contigs"] + " " + work_dir)
    config["contigs"] = os.path.join(config['work_dir'], os.path.basename(config["contigs"]))
    config["sam_file"] = os.path.join(config['work_dir'], "tmp.sam")

    if "bowtie2" in config:
        log.info("bowtie2 found, aligning")
        run_bowtie2(log)
    else:
        log.info("running bwa")
        run_bwa(log)


def run_bowtie2(log):
    global config

    tmp_dir = os.path.join(config["work_dir"], "tmp")
    devnull = "/dev/null"
    universal_sys_call([config["bowtie2"], "-build", config["contigs"], tmp_dir], log, devnull)
    #os.system(config["bowtie2"] + "-build " + config["contigs"] + " " + config["work_dir"] + "tmp > /dev/null")

    universal_sys_call([config["bowtie2"], "-x", tmp_dir, "-1", config["reads1"], "-2",
                     config["reads2"], "-S", config["sam_file"], "-p", str(config["t"]),
                     "--local", "--non-deterministic"], log)
    #os.system(config["bowtie2"] + " -x" + config["work_dir"] + "tmp  -1 " + config["reads1"] + " -2 " + config["reads2"] + " -S " + config["sam_file"] + " -p " + str(config["t"])+ " --local  --non-deterministic")


def run_bwa(log):
    global config

    tmp1_sai_filename = os.path.join(config["work_dir"], "tmp1.sai")
    tmp2_sai_filename = os.path.join(config["work_dir"], "tmp2.sai")
    tmp_sam_filename = os.path.join(config["work_dir"], "tmp.sam")
    isize_txt_filename = os.path.join(config["work_dir"], "isize.txt")
    universal_sys_call([config["bwa"], "index", "-a", "is", config["contigs"], "2"], log)
    universal_sys_call([config["bwa"], "aln", config["contigs"], config["reads1"], "-t",
                       str(config["t"]), "-O", "7", "-E", "2", "-k", "3", "-n", "0.08", "-q", "15"], log, tmp1_sai_filename)
    universal_sys_call([config["bwa"], "aln", config["contigs"], config["reads2"], "-t",
                       str(config["t"]), "-O", "7", "-E", "2", "-k", "3", "-n", "0.08", "-q", "15"], log, tmp2_sai_filename)
    universal_sys_call([config["bwa"], "sampe", config["contigs"], tmp1_sai_filename,
                       tmp2_sai_filename, config["reads1"], config["reads2"]],
                       None, tmp_sam_filename, isize_txt_filename)

    #os.system(config["bwa"] + " index -a is " + config["contigs"] + " 2")
    #os.system(config["bwa"] + " aln  "+ config["contigs"] +" " + config["reads1"] + " -t " + str(config["t"]) + "  -O 7 -E 2 -k 3 -n 0.08 -q 15 >"+config["work_dir"]+ "tmp1.sai" )
    #os.system(config["bwa"] + " aln  "+ config["contigs"] +" " + config["reads2"] + " -t " + str(config["t"]) + "  -O 7 -E 2 -k 3 -n 0.08 -q 15 >"+config["work_dir"]+ "tmp2.sai" )
    #os.system(config["bwa"] + " sampe "+ config["contigs"] +" " + config["work_dir"]+ "/tmp1.sai "+config["work_dir"]+ "/tmp2.sai " + config["reads1"] + " " + config["reads2"] + ">"+config["work_dir"]+ "/tmp.sam 2>"+config["work_dir"]+ "/isize.txt")


def parse_profile(args, log):
    global config

    long_options = "threads= sam-file= output-dir= bwa= contigs= mate-weight= split-dir= bowtie2= 12= insert-size= help debug use-quality use-multiple-aligned skip-masked".split()
    short_options = "1:2:o:s:S:c:t:m:q"

    reads1 = []
    reads2 = []
    reads_mixed = []

    def check_file(f, type, log):
        f = os.path.abspath(f)
        if not os.path.isfile(f):
            error("file with %s (%s) doesn't exist! " % (type, f), log)
        return f

    options, contigs_fpaths = getopt.gnu_getopt(args, short_options, long_options)
    for opt, arg in options:
        if opt in ('-o', "--output-dir"):
            config["output_dirpath"] = os.path.abspath(arg)
            config["make_latest_symlink"] = False
        if opt in ('-c', "--contigs"):
            config["contigs"] = check_file(arg, "contigs", log)
        if opt == '-1':
            reads1.append(check_file(arg, 'left reads', log))
        if opt == '-2':
            reads2.append(check_file(arg, 'right reads', log))
#            config["reads2"] = os.path.abspath(arg)
#            if not os.path.exists(config["reads2"]):
#                log.info("FILE WITH READS DOES NOT EXIST!")
#                usage()
#                sys.exit(1)
        if opt == '--12':
            reads_mixed.append(check_file(arg, 'interleaved left and right reads', log))
#            config["reads_mixed"] = os.path.abspath(arg)
#            if not os.path.exists(config["reads_mixed"]):
#                log.info("FILE WITH READS DOES NOT EXIST!")
#                usage()
#                sys.exit(1)
        if opt == "--bwa":
            config["bwa"] = os.path.abspath(arg)
        if opt in ('-t', "--threads"):
            config["t"] = int(arg)
        if opt in ('-m', "--mate-weight"):
            config["mate_weight"] = float(arg)
        if opt in ('-s', "--sam-file"):
            config["sam_file"] = os.path.abspath(arg)
        if opt in ('-S', "--split-dir"):
            config["split_dir"] = os.path.abspath(arg)
        if opt in ('-q', "--use-quality"):
            config["use_quality"]= 1
        if opt == "--bowtie2":
            if arg != "bowtie2":
                arg = os.path.abspath(arg)
            config["bowtie2"]= arg
        if opt == "--debug":
            config["debug"] = 1
        if opt in ("--use-multiple-aligned"):
            config["use_multiple_aligned"] = 1
        if opt == "--skip-masked":
            config["skip_masked"] = 1
        if opt == "--insert-size":
            config["insert_size"] = int(arg)

    if len(reads1) != len(reads2):
        error("the number of files with left paired reads is not equal to the"
              " number of files with right paired reads!", log)

    if reads1:
        config["reads1"] = reads1
    if reads2:
        config["reads2"] = reads2
    if reads_mixed:
        config["reads_mixed"] = reads_mixed

    work_dir = os.path.join(config["output_dirpath"], "mismatch_corrector_tmp")
    config["work_dir"] = work_dir
    if os.path.isdir(work_dir):
        shutil.rmtree(work_dir)
    os.makedirs(work_dir)
#    os.system ("mkdir -p " + work_dir)
#    os.system("rm -rf " + work_dir + "/*")


def init_config():
    global config
    config = {}

    now = datetime.datetime.now()
    config["output_dirpath"] = "corrector.output." + now.strftime("%Y.%m.%d_%H.%M.%S")+"/"
    config["bwa"] = "bwa"
    config["t"] = int(4)
    config["mate_weight"] = float(1)
    config["use_quality"] = 0
    config["debug"] = 0
    config["use_multiple_aligned"] = 0
    config["skip_masked"] = 0
    config["insert_size"] = int(400)


def process_contig(files):
    log = logging.getLogger('spades')
    samfilename = files[0]
    contig_file = files[1]
    #contigs are always .fasta
    contig_name = contig_file[0:-6]
    logFileName =   contig_name + '.stdout'
    logFile = open(logFileName, 'w')
    ntime = datetime.datetime.now()
    starttime = ntime
    os.system ('./build/release/bin/corrector ' + samfilename + ' ' + contig_file)

    logFile.write(ntime.strftime("%Y.%m.%d_%H.%M.%S") + ": All done. ")
    stime = ntime - starttime
    logFile.write("Time spent: " + str(stime) + "\n")

    logFile.write("Finished processing "+ str(contig_file)) # + ". Used " + str(total_reads) + " reads.\n")
    #logFile.write("replaced: " + str(replaced) + " deleted: "+ str(deleted) +" inserted: " + str(inserted) +'\n')
    logFile.close()

#    return inserted, replaced


def main(args, joblib_path, log=None):

    if len(args) < 1:
        usage()
        sys.exit(0)
    addsitedir(joblib_path)

    os.system ("make -C build/release/corrector/ install=local")
    init_config()
    parse_profile(args, log)

    if not log:
        log = logging.getLogger('spades')
        log.setLevel(logging.DEBUG)

        console = logging.StreamHandler(sys.stdout)
        console.setFormatter(logging.Formatter('%(message)s'))
        console.setLevel(logging.DEBUG)
        log.addHandler(console)

        log_filename = os.path.join(config["output_dirpath"], "corrector.log")
        log_handler = logging.FileHandler(log_filename, mode='w')
        log.addHandler(log_handler)

    log.info("Config: " + str(config))
    if "split_dir" not in config:
#        print "no split dir, looking for sam file"
        if "sam_file" not in config:
            log.info("no sam file, running aligner")
            run_aligner(log)
        else:
            log.info("sam file was found")
            tmp_sam_file_path = os.path.join(config["work_dir"], "tmp.sam")
            shutil.copy2(config["sam_file"], tmp_sam_file_path) # Note: shutil.copy2 is similar to the Unix command cp -p
            #os.system("cp -p "+ config["sam_file"] +" " + config["work_dir"]+"tmp.sam")
            config["sam_file"] = tmp_sam_file_path

        path_to_bin = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../build/release/bin/corrector')
        path_to_config = os.path.join(os.path.dirname(os.path.realpath(__file__)) , '../../configs/corrector/corrector.info')
        print path_to_config
        os.system (path_to_bin + ' ' + config["sam_file"] + ' ' + config["contigs"] + ' ' + ' 0 ' + config["output_dirpath"] + " " + path_to_config)
    #    now = datetime.datetime.now()
    #    res_directory = "corrector.output." + now.strftime("%Y.%m.%d_%H.%M.%S")+"/"


if __name__ == '__main__':
    joblib_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../ext/src/python_libs')
    main(sys.argv[1:], joblib_path)
    
