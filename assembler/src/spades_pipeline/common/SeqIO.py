############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import codecs
import gzip
import sys
import os

fasta_ext = ['.fa', '.fas', '.fasta', '.seq', '.fsa', '.fna', '.ffn', '.frn']
fastq_ext = ['.fq', '.fastq']


def Open(f, mode):
    if f.lower().endswith((".gz", ".gzip")):
        return codecs.getreader('UTF-8')(gzip.open(f, mode))
    else:
        return codecs.open(f, mode, encoding='utf-8')


class Reader:
    def __init__(self, handler):
        self.handler = handler
        self.cash = None

    def FillCash(self):
        if self.cash == None:
            self.cash = self.handler.readline()

    def TrashCash(self):
        self.cash = None

    def Top(self):
        self.FillCash()
        return self.cash

    def readline(self):
        self.FillCash()
        result = self.Top()
        self.TrashCash()
        return result

    def ReadUntill(self, f):
        result = []
        while True:
            next = self.Top()
            if next == "" or f(next):
                return "".join(result)
            self.TrashCash()
            result.append(next.strip())
        return "".join(result)

    def ReadUntillFill(self, buf_size):
        cnt = 0
        result = []
        while not self.EOF() and self.Top() != "" and cnt + len(self.Top().strip()) <= buf_size:
            result.append(self.Top().strip())
            cnt += len(self.Top().strip())
            self.TrashCash()
        if cnt != buf_size:
            raise Exception('The sequence and quality strings for one of reads have different length in file {FILE}. '
                            'Please check the correctness of the file')

        return "".join(result)

    def EOF(self):
        return self.Top() == ""


class SeqRecord:
    def __init__(self, seq, id, qual=None):
        if qual != None and len(qual) != len(seq):
            sys.stdout.write("oppa" + id + "oppa")
        assert qual == None or len(qual) == len(seq)
        self.id = id
        self.seq = seq
        self.qual = qual

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, key):
        return self.seq[key]

    def QualSubseq(self, l, r):
        if self.qual != None:
            return self.qual[l: r]
        return None

    def subseq(self, l, r):
        if l != 0 or r != len(self.seq):
            return SeqRecord(self.seq[l:r], self.id + "(" + str(l + 1) + "-" + str(r) + ")", self.QualSubseq(l, r))
        else:
            return self


def parse(handler, file_type):
    assert file_type in ["fasta", "fastq"]
    if file_type == "fasta":
        return parse_fasta(handler)
    if file_type == "fastq":
        return parse_fastq(handler)


def parse_fasta(handler):
    reader = Reader(handler)
    while not reader.EOF():
        rec_id = reader.readline().strip()
        if len(rec_id) < 1 or rec_id[0] != '>':
            raise Exception("{FILE} is not a valid FASTA file")
        rec_seq = reader.ReadUntill(lambda s: s.startswith(">"))
        yield SeqRecord(rec_seq, rec_id[1:])


def parse_fastq(handler):
    reader = Reader(handler)
    while not reader.EOF():
        rec_id = reader.readline().strip()
        if len(rec_id) < 1 or rec_id[0] != '@':
            raise Exception("{FILE} is not a valid FASTQ file")
        rec_seq = reader.ReadUntill(lambda s: s.startswith("+"))
        tmp = reader.readline()
        if len(tmp) < 1 or (tmp[0] != '+'):
            raise Exception("{FILE} is not a valid FASTQ file")
        rec_qual = reader.ReadUntillFill(len(rec_seq))
        yield SeqRecord(rec_seq, rec_id[1:], rec_qual)


def write(rec, handler, file_type):
    if file_type == "fasta":
        handler.write(">" + rec.id + "\n")
        handler.write(rec.seq + "\n")
    elif file_type == "fastq":
        handler.write("@" + rec.id + "\n")
        handler.write(rec.seq + "\n")
        handler.write("+" + "\n")
        handler.write(rec.qual + "\n")


def FilterContigs(input_handler, output_handler, f, file_type):
    for contig in parse(input_handler, file_type):
        if f(contig):
            write(contig, output_handler, file_type)


def RemoveNs(input_handler, output_handler):
    for contig in parse(input_handler, "fasta"):
        l = 0
        while l < len(contig) and contig[l] == 'N':
            l += 1
        r = len(contig)
        while r > l and contig[r - 1] == 'N':
            r -= 1
        if r > l:
            write(SeqRecord(contig.seq[l:r], contig.id))


def check_extension(fname, extension_list):
    fname = fname.lower()
    # "reads.fastq", ".gz"
    basename_plus_innerext, outer_ext = os.path.splitext(fname)
    if outer_ext not in ['.gz', '.gzip']:
        basename_plus_innerext, outer_ext = fname, ''  # not an archive

    # "reads", ".fastq"
    basename, inner_ext = os.path.splitext(basename_plus_innerext)
    return inner_ext in extension_list


def is_fasta(file_name):
    return check_extension(file_name, fasta_ext)


def is_fastq(file_name):
    return check_extension(file_name, fastq_ext)


def is_bam(file_name):
    return check_extension(file_name, ['.bam'])


def get_read_file_type(file_name):
    if is_fastq(file_name):
        return 'fastq'
    elif is_fasta(file_name):
        return 'fasta'
    elif is_bam(file_name):
        return 'bam'
    else:
        return None
