import bio.core.fasta;
import std.stdio, std.conv, std.string, std.algorithm;

void main(string[] args) {
    if (args.length < 3) { stderr.writeln("Usage: ", args[0], " <reference.fasta> <N>-<M>"); return; }
    auto record = fastaRecords(args[1]).front;
    auto hyphen_pos = std.string.indexOf(args[2], '-');
    auto start = to!size_t(args[2][0 .. hyphen_pos]) - 1;
    auto end = to!size_t(args[2][hyphen_pos + 1 .. $]);
    writeln(">", record.header);
    auto seq = record.sequence[start .. end];
    while (seq.length > 0) {
        writeln(seq[0 .. min(80, $)]);
        seq = seq[min(80, $) .. $];
    }
}
