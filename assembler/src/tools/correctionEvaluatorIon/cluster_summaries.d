import bio.bam.reader;
import std.stdio, std.conv, std.string, std.range;
import std.algorithm : map, filter, sort, uniq, joiner;
import std.array : array;

struct KMer {
    BamRead read;
    uint count;
    double quality;
    alias read this;

    this(BamRead read) {
        this.read = read;
        auto name = read.name[15 .. $];
        count = to!uint(name[0 .. indexOf(name, '-')]);
        quality = to!double(name[lastIndexOf(name, "qual") + 5 .. $]);
    }
}

struct Cluster {
    KMer[] kmers;
    string index;

    void printSummary() {
        auto high_count_kmers = kmers.filter!"a.count > 50"();
        writeln("cluster #", index, " - ", kmers.length, " kmers",
                " (", high_count_kmers.walkLength(), " with count > 50)");
        writeln("    kmers with count > 50 (sequence, count, position, CIGAR, MD tag, edit distance)");
        foreach (kmer; high_count_kmers) {
            writeln("        ", leftJustify(to!string(kmer.sequence), 40), "\t", 
                                kmer.count, "\t",
                                kmer.position, "\t",
                                kmer.cigarString(), "\t", 
                                kmer.read["MD"].is_nothing ? "*" : to!string(kmer.read["MD"]), "\t",
                                kmer.read["NM"].is_nothing ? "*" : to!string(kmer.read["NM"]));
        }
    }
}

void main(string[] args) {
    if (args.length < 2) { stderr.writeln("Usage: ", args[0], " <input.bam>"); return; }
    auto bam = new BamReader(args[1]);

    Cluster cluster;
    foreach (read; bam.reads) {
        auto name = read.name;
        auto index = name[0 .. indexOf(name, '.')];
        if (index == cluster.index)
            cluster.kmers ~= KMer(read);
        else {
            if (cluster.index != "")
                cluster.printSummary();
            cluster.kmers = [KMer(read)];
            cluster.index = index;
        }
    }
}
