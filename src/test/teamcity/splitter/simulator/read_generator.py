#!/usr/bin/env python3
import argparse
import gzip
import math
import os
import random

import gfa_constructor


def get_rc(sequence):
    return "".join([gfa_constructor.RC[nuc] for nuc in sequence[::-1]])


class ReadGenerator(object):
    def __init__(self, read_density, barcode_density, read_length,
                 insert_left, insert_right, barcode_length):
        self.read_density = read_density
        self.barcode_density = barcode_density
        self.read_length = read_length
        self.insert_left = insert_left
        self.insert_right = insert_right
        self.barcode_length = barcode_length
        self.id_seq = "@"
        self.rg_z = "RG:Z:None"

    def _get_fastq_record(self, sequence, id_str):
        return "\n".join([id_str, sequence, "+", "E" * len(sequence)])

    def _get_pair(self, sequence, start_pos, insert_size):
        left_start = start_pos
        left_end = start_pos + self.read_length
        right_start = left_end + insert_size
        right_end = right_start + self.read_length
        left = sequence[left_start:left_end]
        right = get_rc(sequence[right_start:right_end])
        return left, right, left_start, left_end, right_start, right_end

    def _pick_barcode(self, barcodes):
        if random.random() < self.barcode_density:
            return random.choice(barcodes)
        return ""

    def _get_fastq_pair(self, left, right, barcode, path_name,
                        left_start, left_end, right_start, right_end):
        # Both mates share an id carrying path + both coordinate ranges so
        # downstream tooling can re-attribute either read to its source path.
        id_str = "{}{}_L{}-{}_R{}-{} {}".format(
            self.id_seq, path_name,
            left_start, left_end, right_start, right_end,
            self.rg_z)
        if barcode:
            id_str += " BC:Z:" + barcode
        return (self._get_fastq_record(left, id_str),
                self._get_fastq_record(right, id_str))

    def generate_reads(self, sequence, barcodes, path_name):
        assert len(sequence) > 2 * (self.read_length + self.insert_right)
        step = max(1, int(1 / self.read_density))
        # Cap i so right_end = i + 2*read_length + insert_size never exceeds len(sequence).
        stop = len(sequence) - 2 * self.read_length - self.insert_right
        for i in range(0, stop, step):
            insert = random.randint(self.insert_left, self.insert_right)
            left, right, ls, le, rs, re = self._get_pair(sequence, i, insert)
            barcode = self._pick_barcode(barcodes)
            yield self._get_fastq_pair(left, right, barcode, path_name,
                                       ls, le, rs, re)


def parse_fasta(fasta_path):
    sequences = {}
    current_id = None
    current_chunks = []
    with open(fasta_path, "r") as handle:
        for line in handle:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_chunks)
                current_id = line[1:].split()[0]
                current_chunks = []
            else:
                current_chunks.append(line)
        if current_id is not None:
            sequences[current_id] = "".join(current_chunks)
    return sequences


def estimate_reads_per_path(seq_len, step, insert_right, read_length):
    span = seq_len - 2 * read_length - insert_right
    if span <= 0:
        return 0
    return (span + step - 1) // step


def createparser():
    parser = argparse.ArgumentParser(
        description="Simulate paired barcoded reads along graph paths.")
    parser.add_argument("--graph", required=True,
                        help="Graph description in gfa_constructor input format")
    parser.add_argument("--fasta", required=True,
                        help="Fasta with path sequences (gfa_constructor output)")
    parser.add_argument("--output", required=True,
                        help="Output directory for reads_1.fastq / reads_2.fastq")
    parser.add_argument("--read-density", type=float, default=0.2)
    parser.add_argument("--barcode-density", type=float, default=0.99)
    parser.add_argument("--read-length", type=int, default=60)
    parser.add_argument("--insert-left", type=int, default=50)
    parser.add_argument("--insert-right", type=int, default=70)
    parser.add_argument("--barcode-collisions", type=float, default=1.05)
    parser.add_argument("--barcode-length", type=int, default=16)
    parser.add_argument("--seed", type=int, default=None)
    return parser


def main():
    args = createparser().parse_args()
    if args.seed is not None:
        random.seed(args.seed)

    _, graph = gfa_constructor.parse_graph(args.graph)
    fasta_sequences = parse_fasta(args.fasta)

    missing = [p.id for p in graph.paths if p.id not in fasta_sequences]
    if missing:
        raise ValueError("Paths missing from fasta: {}".format(", ".join(missing)))

    path_sequences = [(p.id, fasta_sequences[p.id]) for p in graph.paths]

    step = max(1, int(1 / args.read_density))
    total_reads = sum(estimate_reads_per_path(len(seq), step,
                                              args.insert_right, args.read_length)
                      for _, seq in path_sequences)
    # Pool sized so each barcode is shared by ~barcode_collisions reads on average.
    num_barcodes = max(1, math.ceil(total_reads / args.barcode_collisions))
    barcodes = [gfa_constructor.generate_random_sequence(args.barcode_length)
                for _ in range(num_barcodes)]

    generator = ReadGenerator(
        read_density=args.read_density,
        barcode_density=args.barcode_density,
        read_length=args.read_length,
        insert_left=args.insert_left,
        insert_right=args.insert_right,
        barcode_length=args.barcode_length,
    )

    os.makedirs(args.output, exist_ok=True)
    left_rel = "simulated_1.fastq.gz"
    right_rel = "simulated_2.fastq.gz"
    left_path = os.path.join(args.output, left_rel)
    right_path = os.path.join(args.output, right_rel)
    with gzip.open(left_path, "wt") as left_handle, \
            gzip.open(right_path, "wt") as right_handle:
        for path_name, sequence in path_sequences:
            for left_record, right_record in generator.generate_reads(
                    sequence, barcodes, path_name):
                left_handle.write(left_record + "\n")
                right_handle.write(right_record + "\n")

    yaml_path = os.path.join(args.output, "tenx_dataset.yaml")
    with open(yaml_path, "w") as yaml_handle:
        yaml_handle.write(
            '- "left reads":\n'
            '  - "{}"\n'
            '  "orientation": "fr"\n'
            '  "right reads":\n'
            '  - "{}"\n'
            '  "type": "clouds10x"\n'.format(left_rel, right_rel)
        )


if __name__ == "__main__":
    main()
