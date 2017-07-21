############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2015 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

from __future__ import with_statement
import os
import shlex
import shutil
import stat
import subprocess
import sys
import platform
import re
import gzip
import time
import random
from collections import defaultdict
from genericpath import isdir, exists

from os.path import join

import socket
socket.setdefaulttimeout(20.0)
from urllib2 import urlopen
import xml.etree.ElementTree as ET
import urllib
ncbi_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
main_dir = '/Sid/miheenko/articles/GAGE_C/references'
if platform.system() == 'Darwin':
    sed_cmd = "sed -i '' "
else:
    sed_cmd = 'sed -i '

def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s[0])]


def try_send_request(url):
    attempts = 0
    while attempts < 3:
        try:
            request = urlopen(url, timeout=2000)
            break
        except Exception:
            attempts += 1
            if attempts >= 3:
                print('Cannot established internet connection to download reference genomes! '
                         'Check internet connection or run MetaQUAST with option "--max-ref-number 0".')
            return None
    return request.read()

def download_fastq():
    dirs = os.listdir(main_dir)
    from joblib import Parallel, delayed
    Parallel(n_jobs=8)(delayed(download_one_fastq)(main_dir, dir)
             for dir in dirs)

def sra_id_to_sra_path(sra_id):
    ascp_pref = "" # /sra/sra-instant/reads/ByRun/sra/
    res = "anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/"+ sra_id[:3] +"/" + sra_id[:6] + "/"+ sra_id+ "/" + sra_id + ".sra"
    return res

def download_one_fastq(main_dir, sra_id):
    aspera_string = sra_id_to_sra_path(sra_id)
#    if exists(join(main_dir, sra_id+"_1.fastq")):
    if exists(join(main_dir, sra_id+".sra")):
        print sra_id + " already downloaded"
        return
    os.chdir(main_dir)
#    args = ['fastq-dump', '--split-files', sra_id]
    args = ['ascp', '-i', '/usr/local/etc/asperaweb_id_dsa.openssh', '-k', '1', '-T', '-l', '3000m', aspera_string, '.']
#    print args
    print ('Downloading ' + sra_id)
    subprocess.call(args)
    print ('Downloaded ' + sra_id)
    open(join(main_dir, sra_id+".finished"), 'w')
    return


def download_clade_list(outfilename):
#    ncbi_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    outfile = open(outfilename, 'w')


    response = try_send_request(ncbi_url + 'esearch.fcgi?db=sra&term=Enterobacteriaceae+AND+strategy%3Dwgs+AND+illumina&retmax=105000')
    #response = try_send_request(ncbi_url + 'esearch.fcgi?db=sra&term="bacteria"+'
    #                                       'AND+cluster_public[prop]+"biomol%20dna"+AND+"mate%20pair"&retmax=100000')
    xml_tree = ET.fromstring(response)
    if xml_tree.find('Count').text == '0':  # Organism is not found
        return None
    print "Found %s datasets" % xml_tree.find('Count').text
    ref_fpath = None
    ref_ids = [ref.text for ref in xml_tree.find('IdList').findall('Id')]
    sra_by_refs = defaultdict(dict)
    len_by_refs = defaultdict(int)
    for i, assembly_id in enumerate(ref_ids):
        sample_name = None
        request = ncbi_url + 'elink.fcgi?dbfrom=sra&db=taxonomy&id=%s' % assembly_id
        response = try_send_request(request)
        if response is None:
             continue
        xml_tree = ET.fromstring(response)
        if xml_tree is None or not xml_tree:
            print "NoneRequest"
            continue
        if xml_tree.find('LinkSet') is None:
            continue
        if xml_tree.find('LinkSet').find('LinkSetDb') is None:
            continue
        if xml_tree.find('LinkSet').find('LinkSetDb').find('Link') is None:
            continue
        if xml_tree.find('LinkSet').find('LinkSetDb').find('Link').find('Id') is None:
            continue

        tax_id = xml_tree.find('LinkSet').find('LinkSetDb').find('Link').find('Id').text
        outfile.write(assembly_id+" "+ tax_id + '\n')

    return


def download_ref(ref_name, assembly_id):
    ncbi_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    sample_dir = join(main_dir, ref_name)
    if not isdir(sample_dir):
        os.mkdir(sample_dir)
    ref_fpath = join(sample_dir, ref_name.replace(' ', '_') + '.fasta')
    if exists(ref_fpath) and os.path.getsize(ref_fpath) > 100:
        print (ref_fpath + ' already exists')
        return ref_fpath
    response = try_send_request(
        ncbi_url + 'elink.fcgi?dbfrom=assembly&db=nuccore&linkname="assembly_nuccore_refseq"&id=%s' % assembly_id)
    xml_tree = ET.fromstring(response)
    if xml_tree is None:
        return None

    link_set = xml_tree.find('LinkSet')
    if link_set is None or xml_tree.find('LinkSet').find('LinkSetDb') is None:
        response = try_send_request(
            ncbi_url + 'elink.fcgi?dbfrom=assembly&db=nuccore&linkname="assembly_nuccore_insdc"&id=%s' % assembly_id) ## try in genbank
            
        xml_tree = ET.fromstring(response)
        if xml_tree is None:
            return None

        link_set = xml_tree.find('LinkSet')
        if link_set is None or xml_tree.find('LinkSet').find('LinkSetDb') is None:
            return None

    is_first_piece = False
    fasta_files = []
    ref_ids = [link.find('Id').text for link in xml_tree.find('LinkSet').find('LinkSetDb').findall('Link')]
    for ref_id in ref_ids:
        fasta = try_send_request(ncbi_url + 'efetch.fcgi?db=sequences&id=%s&rettype=fasta&retmode=text' % ref_id)
        if fasta:
            fasta_files.append(fasta)
    fasta_names = [f.split('|')[-1] for f in fasta_files]
    k = 0
    with open(ref_fpath, "w") as fasta_file:
        for name, fasta in zip(fasta_names, fasta_files):
            if not is_first_piece:
                is_first_piece = True
            else:
                fasta = '\n' + fasta.rstrip()
            fasta_file.write(fasta.rstrip())
            k += 1
            print (ref_fpath + ' is downloaded ' + str(k) + ' from ' + str(len(fasta_files)))
    return ref_fpath


def show_progress(a, b, c):
    if a > 0 and a % int(c/(b*100)) == 0:
        print("% 3.1f%% of %d bytes\r" % (min(100, int(float(a * b) / c * 100)), c)),
        sys.stdout.flush()

def check_sra_id(sra_id):
    sample_name = None
    response = try_send_request(ncbi_url + 'efetch.fcgi?db=sra&id=%s&tool=quast&email=quast.support@bioinf.spbau.ru' % sra_id)
    if not response:
        print "no response"
        return False
    xml_tree = ET.fromstring(response)
    if xml_tree.find('EXPERIMENT_PACKAGE') is None or xml_tree.find('EXPERIMENT_PACKAGE').find('SAMPLE') is None:
        return False
    sample_id = xml_tree.find('EXPERIMENT_PACKAGE').find('SAMPLE').find('SAMPLE_NAME').find('TAXON_ID').text
    sample_name_el = xml_tree.find('EXPERIMENT_PACKAGE').find('SAMPLE').find('SAMPLE_NAME').find('SCIENTIFIC_NAME')
    if sample_name_el is None:
        print 'no info for ' + sra_id
        return False
    if xml_tree.find('EXPERIMENT_PACKAGE').find('EXPERIMENT').find('PLATFORM').find('ILLUMINA') is None:
        print 'no illumina reads for ' + sra_id
        return False
    sample_name = sample_name_el.text
    sample_name = re.sub('[/.= ]', '_', sample_name)
    title = xml_tree.find('EXPERIMENT_PACKAGE').find('STUDY').find('DESCRIPTOR').find('STUDY_TITLE').text
    library = xml_tree.find('EXPERIMENT_PACKAGE').find('EXPERIMENT').find('DESIGN').find('LIBRARY_DESCRIPTOR').find('LIBRARY_LAYOUT').getchildren()[0]
    library_len = '0'
    if 'NOMINAL_LENGTH' in library.attrib:
        library_len = library.attrib['NOMINAL_LENGTH']
    if xml_tree.find('EXPERIMENT_PACKAGE') is None or xml_tree.find('EXPERIMENT_PACKAGE').find('RUN_SET') is None or xml_tree.find('EXPERIMENT_PACKAGE').find('RUN_SET').findall('RUN') is None:
        print "empty tree in trange  place for" + sra_id
        return False
    runs = xml_tree.find('EXPERIMENT_PACKAGE').find('RUN_SET').findall('RUN')
    for run in runs:
        try:
            size = run.attrib['size'] if 'size' in run.attrib else 0
        except:
            size = 0
        strain_name = ''
        info = ''
        tmp = run.find('IDENTIFIERS')
        if tmp:
            tmp2 = run.find('IDENTIFIERS').find('SUBMITTER_ID')
            if tmp2:
                info = run.find('IDENTIFIERS').find('SUBMITTER_ID').text
        if info:
            strain_name = re.findall(r'strain (\S+)', info)
            if strain_name:
                sample_name += '_strain_' + strain_name[0]
        if int(size) / 1024 / 1024 < 100:  # reads is smaller than 400 Mb
            print sample_name_el.text + ' too small sra ' + sra_id
            return False
        sra_found = True
        sra_id2 = run.find('IDENTIFIERS').find('PRIMARY_ID').text
        #sra_by_refs[sample_name][library.tag] = (size, sra_id2, library_len, title)
        #sample_name = None
        print sra_id
        print sample_name, library.tag, size, sra_id2
    return True

def parse_list(listname, sra_id_list):
#    ncbi_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    random.seed(239)
    file = open(listname, 'r')
    sra_file = open(sra_id_list, 'w')
    arr = {}
    for line in file:
        sample_id, tax_id = line.strip().split()
        if not check_sra_id(sample_id):
            continue
        if tax_id in arr.keys():
            arr[tax_id].append(sample_id)
        else:
            arr[tax_id] = [sample_id]
    for tax in sorted(arr):
        print "%s: %s " %(tax, len(arr[tax]))
    for tax in arr:
        tmp_id =  random.randint(0, len(arr[tax]) -1 )
        print tmp_id, len(arr[tax])
        sample_sra_id = arr[tax][tmp_id]
        sra_file.write(sample_sra_id + " " + tax + "\n")
    return


def download_sra_from_list(sra_id_list, download_dir):
#    download_dir = '/Bmo/dantipov/jgi_5625_7625/'

    sra_file = open(sra_id_list, 'r')
    sra_by_refs = defaultdict(dict)
    for line in sra_file:
#        sra_id, tax_id = line.strip().split()
        sra_id = line.strip()
        sample_name = None
        response = try_send_request(ncbi_url + 'efetch.fcgi?db=sra&id=%s&tool=quast&email=quast.support@bioinf.spbau.ru' % sra_id)
        if not response:
            print 'no responce'
            continue
        xml_tree = ET.fromstring(response)
        if xml_tree.find('EXPERIMENT_PACKAGE') is None or xml_tree.find('EXPERIMENT_PACKAGE').find('SAMPLE') is None:
            print 'no experimient package r no sample'
            continue
        sample_id = xml_tree.find('EXPERIMENT_PACKAGE').find('SAMPLE').find('SAMPLE_NAME').find('TAXON_ID').text
        sample_name_el = xml_tree.find('EXPERIMENT_PACKAGE').find('SAMPLE').find('SAMPLE_NAME').find('SCIENTIFIC_NAME')
        if sample_name_el is None:
            print 'no info for ' + sra_id
            continue
        if xml_tree.find('EXPERIMENT_PACKAGE').find('EXPERIMENT').find('PLATFORM').find('ILLUMINA') is None:
            print 'no illumina reads for ' + sra_id
            continue
        sample_name = sample_name_el.text
        sample_name = re.sub('[/.= ]', '_', sample_name)
        title = xml_tree.find('EXPERIMENT_PACKAGE').find('STUDY').find('DESCRIPTOR').find('STUDY_TITLE').text
        library = xml_tree.find('EXPERIMENT_PACKAGE').find('EXPERIMENT').find('DESIGN').find('LIBRARY_DESCRIPTOR').find('LIBRARY_LAYOUT').getchildren()[0]
        library_len = '0'
        if 'NOMINAL_LENGTH' in library.attrib:
            library_len = library.attrib['NOMINAL_LENGTH']
        runs = xml_tree.find('EXPERIMENT_PACKAGE').find('RUN_SET').findall('RUN')
        for run in runs:
            try:
                size = run.attrib['size'] if 'size' in run.attrib else 0
            except:
                size = 0
            strain_name = ''
            info = ''
            tmp = run.find('IDENTIFIERS')
            if tmp:
                tmp2 = run.find('IDENTIFIERS').find('SUBMITTER_ID')
                if tmp2:
                    info = run.find('IDENTIFIERS').find('SUBMITTER_ID').text
            if info:
                strain_name = re.findall(r'strain (\S+)', info)
                if strain_name:
                    sample_name += '_strain_' + strain_name[0]
            if int(size) / 1024 / 1024 < 100:  # reads is smaller than 100 Mb
                print sample_name_el.text.encode('utf-8') + ' too small sra ' + sra_id
                continue
            sra_found = True
            sra_id2 = run.find('IDENTIFIERS').find('PRIMARY_ID').text
            sra_by_refs[sample_name][library.tag] = (size, sra_id2, library_len, title)
            #sample_name = None
#            print sra_id, tax_id            
            print sample_name.encode('utf-8'), library.tag, size, sra_id2, library_len, title.encode('utf-8')
            download_one_fastq(download_dir, sra_id2)
    sra_file.close()
            
def dump_one_fastq(l, dir):
    arr = l.strip().split('.')
    if len(arr) != 2:
        return
    ext = l.strip().split('.')[-1]
    filename = l.strip().split('.')[0]
    if ext != "sra":
        print ext + " is wrong ext"
        return
    if exists(os.path.join(dir,filename + "_1.fastq")):
        print filename + "already dumped"
    else:
        args = ['fastq-dump', '--split-3', os.path.join(dir, l), '-O', dir]
        print args
        subprocess.call(args)
        print ('Dumped ' + l)

def dump_all_fastq(dir):
    from joblib import Parallel, delayed
    Parallel(n_jobs=16)(delayed(dump_one_fastq)(l, dir)
             for l in os.listdir(dir))

def main():
    #logger.print_timestamp()
#    tax_to_sra_list= sys.argv[1]
    sra_id_list = sys.argv[1]
#    download_clade_list(tax_to_sra_list)
    #parse_list(tax_to_sra_list, sra_id_list)
    work_dir = sys.argv[2]
#    dump_all_fastq(work_dir)
    for i in range (1, 100):
        print "downloading restarted"
        download_sra_from_list(sra_id_list, work_dir)
        print "downloading attempt ended"
        dump_all_fastq(work_dir)
        print "dumping ended"
        time.sleep(1100)    
        
    #download_fastq()

if __name__ == '__main__':
    main()


