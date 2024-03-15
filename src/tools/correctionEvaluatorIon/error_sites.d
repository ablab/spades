import bio.bam.reader, bio.bam.baseinfo;

import std.stdio, std.algorithm, std.range, std.traits, std.array;

struct Bin {
    size_t[] positions;
}

Bin[] binned(T)(T numbers, int max_adj_diff) {
    if (numbers.empty) return [];
    Bin[] result;
    Bin current_bin;
    current_bin.positions = [numbers[0]];
    foreach (num; numbers[1 .. $]) {
        if (num <= current_bin.positions.back + max_adj_diff) {
            current_bin.positions ~= num;
        } else {
            result ~= current_bin;
            current_bin = Bin.init;
            current_bin.positions = [num];
        }
    }
    if (current_bin.positions.length > 0)
        result ~= current_bin;
    return result;
}

void main(string[] args) {
    if (args.length < 2) { stderr.writeln("Usage: ", args[0], " <input.bam>"); return; }
    auto bam = new BamReader(args[1]);

    uint[size_t] error_positions;    

    foreach (read; bam.reads) {
        alias basesWith!("MD", Option.cigarExtra) baseInfo;
        ElementType!(ReturnType!(baseInfo!BamRead))[2048] bases_buf;
        size_t sz;
        foreach (base; baseInfo(read))
            bases_buf[sz++] = base;
        auto bases = bases_buf[0 .. sz];

        foreach (baseinfo; bases) {
            with (baseinfo) {
                bool error = cigar_operation.type == 'M' && base != reference_base;

                if (read.strand == '+') {
                    error |= cigar_operation.type == 'I' && cigar_operation_offset == 0;
                    error |= cigar_operation_offset == cigar_operation.length - 1 
                             && !cigar_after.empty && cigar_after.front.type == 'D';
                } else {
                    error |= cigar_operation.type == 'I' && 
                             cigar_operation_offset == cigar_operation.length - 1;
                    error |= cigar_operation_offset == 0 && !cigar_before.empty
                             && cigar_before.back.type == 'D';
                }

                if (error)
                    error_positions[position] += 1;
            }
        }
    }

    auto bins = error_positions.keys.sort.binned(10);

    foreach (bin; bins)
    {
        size_t n_reads;
        foreach (pos; bin.positions)
            n_reads += error_positions[pos];
        if (bin.positions.length == 1)
            writeln(bin.positions.front, "\t", n_reads);
        else
            writeln(bin.positions.front, "-", bin.positions.back, "\t", n_reads);
    }
}
