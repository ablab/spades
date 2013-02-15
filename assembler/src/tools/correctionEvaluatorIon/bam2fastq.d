import bio.bam.reader, bio.core.sequence;
import std.stdio, std.range, std.algorithm, std.array;

bool is_not_desirable(char cigar_op_type) @property {
  return cigar_op_type == 'S' || cigar_op_type == 'I';
}

void main(string[] args) {
  auto bam = new BamReader(args[1]);
  auto f = File(bam.filename ~ ".fastq", "w+");

  foreach (read; bam.reads().filter!"!a.is_unmapped"()) {
    with (read) {
      auto seq = nucleotideSequence(sequence);
      auto qual = base_qualities.map!"cast(char)(a + 33)"();

      auto cigar = read.cigar;
      if (cigar.empty)
        continue;

      // get rid of hard clips if they exist
      if (cigar.front.type == 'H')
        cigar.popFront();
      if (cigar.back.type == 'H')
        cigar.popBack();

      // remove soft clips and insertions from read ends
      // so as to make comparison easier
      while (!cigar.empty && cigar.front.type.is_not_desirable) {
        seq.popFrontN(cigar.front.length);
        qual.popFrontN(cigar.front.length);
        cigar.popFront();
      }

      while (!cigar.empty && cigar.back.type.is_not_desirable) {
        seq.popBackN(cigar.back.length);
        qual.popBackN(cigar.back.length);
        cigar.popBack();
      }

      // don't output read if it has very long indel in the middle --
      // this makes it hard to compare it with the corrected one
      bool good = true;
      foreach (op; cigar) {
        if (!op.is_match_or_mismatch && op.length > 5) {
          good = false;
          break;
        }
      }

      if (!good)
        continue;

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
