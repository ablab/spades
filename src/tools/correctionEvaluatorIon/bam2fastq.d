import bio.bam.reader, bio.core.sequence;
import std.stdio, std.range, std.algorithm, std.array;

void main(string[] args) {
  if (args.length < 2) { stderr.writeln("Usage: ", args[0], " <input.bam>"); return; }
  auto bam = new BamReader(args[1]);
  auto f = File(bam.filename ~ ".fastq", "w+");

  foreach (read; bam.reads().filter!"!a.is_unmapped"()) {
    with (read) {
      auto seq = nucleotideSequence(sequence);
      auto qual = base_qualities.map!"cast(char)(a + 33)"();

      f.write('@', name);
      if (is_paired)
        f.writeln('/', is_first_of_pair ? '1' : '2');
      else
        f.writeln();

      if (strand == '+')
        f.writeln(seq, "\n+\n", qual);
      else
        f.writeln(seq.reverse, "\n+\n", qual.retro());
    }
  }
}
