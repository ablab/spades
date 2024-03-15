import bio.bam.reader;
import bio.bam.read;
import bio.bam.baseinfo;
import bio.core.base;

import std.array;
import std.conv;
import std.traits;
import std.math;
import std.range;
import std.algorithm;
import std.parallelism;
import std.stdio;

immutable K = 5;

class HintTable {

    private {
        uint[][][K + 1] insertions;
        uint[][][K + 1] deletions;
        uint[][][K + 1] mismatches;
        uint[][][K + 1] total;
    }

    this() {
        foreach (k; 1 .. K + 1) {
            auto n_kmers = 1 << (2 * k);
            insertions[k] = new uint[][](n_kmers, n_kmers);
            deletions[k] = new uint[][](n_kmers, n_kmers);
            mismatches[k] = new uint[][](n_kmers, n_kmers);
            total[k] = new uint[][](n_kmers, n_kmers);
        }
    }

    void mergeWith(HintTable other) {
        foreach (k; 1 .. K + 1) {
            auto n_kmers = 1 << (2 * k);
            foreach (kmer1; 0 .. n_kmers) {
                insertions[k][kmer1][] += other.insertions[k][kmer1][];
                deletions[k][kmer1][] += other.deletions[k][kmer1][];
                mismatches[k][kmer1][] += other.mismatches[k][kmer1][];
                total[k][kmer1][] += other.total[k][kmer1][];
            }
        }
    }

    void printTable() {
        writeln("ref\tread\ttotal\tI\tD\tX");
        foreach (k1; 0 .. (1U << (2 * K))) {
            foreach (k2; 0 .. (1U << (2 * K))) {
                auto c = total[K][k1][k2];
                if (c > 0) {
                    writeln(uintToKmer(k1, K), " ", uintToKmer(k2, K), " ", c,
                            "\t", insertions[K][k1][k2],
                            "\t", deletions[K][k1][k2],
                            "\t", mismatches[K][k1][k2]);
                }
            }
        }
    }

    void saveTo(string filename) {
        auto f = File(filename, "w+");
        foreach (k; 1 .. K + 1) {
            auto n_kmers = 1U << (2 * k);
            foreach (i; 0 .. n_kmers) {
                f.rawWrite(insertions[k][i]);
                f.rawWrite(deletions[k][i]);
                f.rawWrite(mismatches[k][i]);
                f.rawWrite(total[k][i]);
            }
        }
    }

    void loadFrom(string filename) {
        auto f = File(filename);
        foreach (k; 1 .. K + 1) {
            auto n_kmers = 1U << (2 * k);
            foreach (i; 0.. n_kmers) {
                f.rawRead(insertions[k][i][]);
                f.rawRead(deletions[k][i][]);
                f.rawRead(mismatches[k][i][]);
                f.rawRead(total[k][i][]);
            }
        }
    }

    void appendFrom(string filename) {
        auto other = new HintTable();
        other.loadFrom(filename);
        this.mergeWith(other);
    }

    enum Hint {
        mismatch,
        insertion,
        deletion
    }

    Hint getHint(size_t k, uint kmer1, uint kmer2) {
        while (k > 3 && total[k][kmer1][kmer2] == 0) {
            k--;
            kmer1 &= (1 << (2 * k)) - 1;
            kmer2 &= (1 << (2 * k)) - 1;
        }
        auto m = mismatches[k][kmer1][kmer2];
        auto i = insertions[k][kmer1][kmer2];
        auto d = deletions[k][kmer1][kmer2];
        if (m > i && m > d)
            return Hint.mismatch;
        if (i >= m && i >= d)
            return Hint.insertion;
        if (d >= m && d >= i)
            return Hint.deletion;
        return Hint.mismatch;
    }

    void saveHintsInc(string filename) {
        auto f = File(filename, "w+");
        uint[][K + 1] tables;
        foreach (k; 1 .. K + 1) {
            auto n_kmers = 1 << (2 * k);
            tables[k] = new uint[(n_kmers * n_kmers) / 16];
            foreach (k1; 0 .. n_kmers) {
                foreach (k2; 0 .. n_kmers) {
                    auto hint = cast(ubyte)getHint(k, k1, k2);
                    auto id = k1 + k2 * n_kmers;
                    auto shift = 2 * (id % 16);
                    tables[k][id / 16] |= (hint << shift);
                }
            }
        }
        f.writefln("%(%sU,%)U", tables[1 .. K + 1].joiner());
    }

    void collectErrorCounts(BamReader bam) {
        foreach (read; bam.reads) {
            if (read.is_unmapped || read["MD"].is_nothing || read.cigar.empty)
                continue;

            auto edit_distance = read["NM"].to!int();
            if (edit_distance == 0)
                continue;

            collectErrorCounts(read);
        }
    }

    private void collectErrorCounts(BamRead read) {
        alias basesWith!"MD" baseInfo;
        ElementType!(ReturnType!(baseInfo!BamRead))[2048] baseinfo_buf = void;
        size_t len = 0;
        foreach (info; baseInfo(read))
            baseinfo_buf[len++] = info;
        auto baseinfo = baseinfo_buf[0 .. len];

        size_t pos = 0;
        while (pos < len && baseinfo[pos].cigar_operation.type == 'S')
            ++pos;

        Base5[K] ref_kmer = void;
        Base5[K] read_kmer = void;

        while (true) {
            size_t ref_kmer_len;
            size_t read_kmer_len;
            while (pos < len && baseinfo[pos].cigar_operation.type == 'M'
                             && baseinfo[pos].base == baseinfo[pos].reference_base)
            {
                pos++;
            }

            if (pos == len || baseinfo[pos].cigar_operation.type == 'S')
                break;

            auto op = baseinfo[pos].cigar_operation;

            for (size_t i = pos; i < len; i++) {
                if (read_kmer_len == K && ref_kmer_len == K) break;

                if (read_kmer_len < K && baseinfo[i].cigar_operation.is_query_consuming)
                    read_kmer[read_kmer_len++] = Base5(baseinfo[i].base);
               
                if (ref_kmer_len < K && baseinfo[i].cigar_operation.is_reference_consuming)
                    ref_kmer[ref_kmer_len++] = Base5(baseinfo[i].reference_base);
            }

            if (read_kmer_len == K && ref_kmer_len == K) {
                auto fst = kmerToUint(ref_kmer);
                auto snd = kmerToUint(read_kmer);
                if (op.type == 'M') {
                    for (size_t k = K; k >= 1; k--) {
                        fst &= (1 << (2 * k)) - 1;
                        snd &= (1 << (2 * k)) - 1;
                        mismatches[k][fst][snd] += 1;
                        mismatches[k][snd][fst] += 1;
                        total[k][fst][snd] += 1;
                        total[k][snd][fst] += 1;
                    }
                } else if (op.type == 'I') {
                    for (size_t k = K; k >= 1; k--) {
                        fst &= (1 << (2 * k)) - 1;
                        snd &= (1 << (2 * k)) - 1;
                        insertions[k][fst][snd] += 1;
                        deletions[k][snd][fst] += 1;
                        total[k][fst][snd] += 1;
                        total[k][snd][fst] += 1;
                    }
                }
            }

            ++pos;

            while (pos < len && !baseinfo[pos].cigar_operation.is_reference_consuming)
                ++pos;
        }
    }

    private static uint kmerToUint(ref Base5[K] kmer) {
        uint num = 0;
        uint i = 0;
        foreach (base; kmer[])
            num |= (base.internal_code << (2 * i++));
        return num;
    }

    private static string uintToKmer(uint kmer, size_t len) {
        string s;
        for (size_t i = 0; i < len; i++)
            s ~=  Base5.fromInternalCode((kmer >> (2 * i)) & 0x3);
        return s;
    }
}

void printUsage() {
    stderr.writeln(`usage: ./ion_dist_table collect <input.bam> <output.ec>`);
    stderr.writeln(`       ./ion_dist_table merge <output.ec> <input1.ec> <input2.ec> ...`);
    stderr.writeln(`       ./ion_dist_table generate <output.inc> <input1.ec> <input2.ec> ...`);
    stderr.writeln(`       ./ion_dist_table view <input.ec>`);
}

void main(string[] args) {
    if (args.length == 1) {
        printUsage();
        return;
    }

    auto hints = new HintTable();
    switch (args[1]) {
        case "collect":
            if (args.length < 4) { printUsage(); return; }
            auto pool = new TaskPool(4);
            scope(exit) pool.finish();
            auto bam = new BamReader(args[2], pool);
            hints.collectErrorCounts(bam);
            hints.saveTo(args[3]);
            return;
        case "merge":
            if (args.length < 4) { printUsage(); return; }
            foreach (fn; args[3 .. $])
                hints.appendFrom(fn);
            hints.saveTo(args[2]);
            return;
        case "generate":
            if (args.length < 4) { printUsage(); return; }
            foreach (fn; args[3 .. $])
                hints.appendFrom(fn);
            hints.saveHintsInc(args[2]);
            return;
        case "view":
            if (args.length < 3) { printUsage(); return; }
            hints.loadFrom(args[2]);
            hints.printTable();
            return;
        default:
            printUsage();
            return;
    }
}
