#!/usr/bin/env python

############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Convert binary output of spades-kmercount to fasta format')
    parser.add_argument("-k",
                        type=int,
                        default=21,
                        help="K-mer length",
                        action="store")
    parser.add_argument("-o",
                        dest="output",
                        default="final_kmers.fasta",
                        help="output file to use",
                        action="store")
    parser.add_argument(dest="input_filename",
                        metavar="filename",
                        help="path to binary file with kmers(spades-kmercount output)",
                        action="store")
    return parser.parse_args()


def main():
    args = parse_args()
    nucl = "ACGT"
    k_len = args.k
    buckets_per_kmer = (k_len*2 - 1)//64 + 1
    bytes_per_bucket = 8
    kmer_cnt = 0
    with open(args.output, 'w') as fw:
        with open(args.input_filename, 'br') as f:
            data = f.read(bytes_per_bucket * buckets_per_kmer)
            while data:
                kmer = ""
                for bt in data:
                    for i in range(4):
                        kmer += nucl[(bt >> (i*2))&3]
                    kmer = kmer[:k_len]
                fw.write(">kmer" + str(kmer_cnt) + "\n")
                fw.write(kmer + "\n")
                kmer_cnt += 1
                data = f.read(bytes_per_bucket * buckets_per_kmer)


if __name__ == "__main__":
    main()