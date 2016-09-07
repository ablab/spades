#!/usr/bin/env python3

import os
import argparse

def parse_args():
	parser = argparse.ArgumentParser(description="Kmers mpl filter")
	parser.add_argument("-om", "--one-min", default=3, type=int, help="min kmer mpl in one sample")
	parser.add_argument("-am", "--all-min", default=3, type=int, help="min kmer mpl in all samples")
	parser.add_argument("-kl", "--kmer-len", default=31, type=int, help="kmer length")
	parser.add_argument("samples_dir", help="directory with samples")
	parser.add_argument("output", help="output files prefix")
	args = parser.parse_args()
	return args

def calc_mpl(args):
	if not os.path.exists(args.samples_dir):
		os.makedirs(args.samples_dir)

	files = [f for f in os.listdir(args.samples_dir) if os.path.isfile(os.path.join(args.samples_dir, f))]

	cmd = "/home/toxa31/work/algorithmic-biology/assembler/src/kmer_count_filter/kmer_count_filter -kl {} -one-min {} -all-min {}".format(
		args.kmer_len, args.one_min, args.all_min)

	for f in files:
		cmd = cmd + " " + args.samples_dir + "/" + f

	cmd = cmd + " " + args.output

	print(cmd)

	os.system(cmd)

def main():
	args = parse_args()
	calc_mpl(args)

main()