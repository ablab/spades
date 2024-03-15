import bio.bam.reader, bio.bam.read, bio.bam.md.reconstruct, bio.bam.baseinfo, bio.bam.writer;
import bio.core.sequence, bio.core.fasta, bio.core.base;
import std.stdio, std.range, std.conv, std.algorithm, 
       std.parallelism, std.array, std.traits, std.exception;

auto mapped_sequence(ref BamRead read) @property {
  auto sequence = read.sequence;
  auto cigar = read.cigar;
  while (!cigar.empty && !cigar.front.is_reference_consuming)
  {
    sequence.popFrontN(cigar.front.length);
    cigar.popFront();
  }

  while (!cigar.empty && !cigar.back.is_reference_consuming)
  {
    sequence.popBackN(cigar.back.length);
    cigar.popBack();
  }
  
  return sequence;
}

// FIXME: more fine-grained statistics about 
//        corrected/miscorrected/uncorrected errors

class Comparator
{
  private {
    BamReader raw_;
    BamReader corrected_;

    BamWriter ruined_;
    BamWriter damaged_;
  }

  this(BamReader raw, BamReader corrected) {
    raw_ = raw;
    corrected_ = corrected;
    
    worsened_ = File("worsened.txt", "w+");
  }

  /// Collect statistics
  void run() {

    ruined_ = new BamWriter("ruined.bam");
    ruined_.writeSamHeader(corrected_.header);
    ruined_.writeReferenceSequenceInfo(corrected_.reference_sequences);
    scope (exit) ruined_.finish();

    damaged_ = new BamWriter("damaged.bam");
    damaged_.writeSamHeader(corrected_.header);
    damaged_.writeReferenceSequenceInfo(corrected_.reference_sequences);
    scope (exit) damaged_.finish();

    auto last_position = reduce!max(0, 
                            raw_.reads().filter!"!a.is_unmapped"
                                .map!"a.position + a.basesCovered() - 1");

    auto raw_reads = raw_.reads;
    foreach (r2; corrected_.reads) {
      if (r2.position >= last_position)
          continue;
      while (!raw_reads.empty && raw_reads.front.name != r2.name)
        raw_reads.popFront();
      assert(!raw_reads.empty);
      auto r1 = raw_reads.front;
      updateStats(r1, r2);
    }
  }

  private void outputDistribution(D)(string filename, auto ref D distribution) {
    auto f = File(filename, "w+");
    foreach (size_t lev_dist, count; distribution)
      if (count > 0)
        f.writeln(lev_dist, '\t', count);
  }

  private static bool isTrimmed(char op_t) {
    return op_t == 'S' || op_t == 'I';
  }

  // Get edit path for aligned read.
  private static EditOp[] cigarPath(BamRead read) {
    EditOp[] path;

    alias basesWith!("MD", Option.cigarExtra) baseInfo;
    ElementType!(ReturnType!(baseInfo!BamRead))[2048] bases_buf;
    size_t sz;
    foreach (base; baseInfo(read))
      bases_buf[sz++] = base;
    auto bases = bases_buf[0 .. sz];

    // it's assumed that insertions/soft clips on ends are trimmed by bam2fastq.d
    while (!bases.empty && isTrimmed(bases.front.cigar_operation.type))
      bases.popFront();
    while (!bases.empty && isTrimmed(bases.back.cigar_operation.type))
      bases.popBack();

    foreach (baseinfo; bases) {
      if (baseinfo.cigar_operation.type == 'M')
        path ~= baseinfo.base == baseinfo.reference_base ? EditOp.none : EditOp.substitute;
      else if (baseinfo.cigar_operation.is_query_consuming)
        path ~= EditOp.insert;
      else
        assert(0);

      if (baseinfo.cigar_operation_offset == baseinfo.cigar_operation.length - 1)
        foreach (op; baseinfo.cigar_after)
          if (op.is_reference_consuming && !op.is_query_consuming)
            foreach (i; 0 .. op.length)
              path ~= EditOp.remove;
          else
            break; // stop at insertion/soft clip
    }

    return path;
  }

  // format read sequence so that it's aligned with edit path
  private static string sequenceStr(BamRead read) {
    auto cigar = read.cigar;
    auto seq = read.sequence;
    while (!cigar.empty && isTrimmed(cigar.front.type)) {
      seq.popFrontN(cigar.front.length);
      cigar.popFront();
    }
    while (!cigar.empty && isTrimmed(cigar.back.type)) {
      seq.popBackN(cigar.back.length);
      cigar.popBack();
    }
  
    char[] res;
    foreach (op; cigar) {
      if (op.is_query_consuming) {
        foreach (i; 0 .. op.length) {
          res ~= seq.front;
          seq.popFront();
        }
      } else if (op.is_reference_consuming) {
        res ~= std.array.replicate("*", op.length);
      }
    }
    if (read.strand == '-') {
      res.reverse();
      foreach (ref c; res)
        if (c != '*')
          c = Base5(c).complement.asCharacter;
    }
    return assumeUnique(res);
  }

  private void updateStats(BamRead r1, BamRead r2) {
    if (r1.is_unmapped || r2.is_unmapped)
      return;

    if (r1.strand != r2.strand)
      return; // comparison is meaningless

    if (r1.strand == '+' && r1.position != r2.position)
        return;

    if (r1.strand == '-' && r1.position + r1.basesCovered() != r2.position + r2.basesCovered())
        return; // one of reads was likely misaligned

    if (r1["NM"].is_nothing || r2["NM"].is_nothing)
      return;

    auto old_cigar_path = cigarPath(r1);
    auto old_lev_dist = to!long(r1["NM"]);
    auto old_bc = r1.basesCovered();

    auto new_cigar_path = cigarPath(r2);
    auto new_lev_dist = to!long(r2["NM"]);
    auto new_bc = r2.basesCovered();

    delta_cov_ += to!long(new_bc) - to!long(old_bc);

    lev_distance_distribution_[new_lev_dist] += 1;

    ++total_;
    delta_ += new_lev_dist - old_lev_dist;
    old_delta_ -= old_lev_dist;

    if (old_lev_dist == 0)
      ++n_perfect_;
    else
      ++n_err_;

    if (old_lev_dist == 1)
      ++n_one_err_;

    bool degraded = old_lev_dist < new_lev_dist;

    if (old_lev_dist == 0) {
      if (new_lev_dist > 0)
      {
        auto r1_seq = r1.mapped_sequence;
        auto r2_seq = r2.mapped_sequence;
        if (r2_seq.length >= r1_seq.length &&
            ((r1.strand == '+' && equal(r1_seq, r2_seq.take(r1_seq.length))) ||
             (r1.strand == '-' && equal(r1_seq.retro(), 
                                        r2_seq.retro().take(r1_seq.length)))))
        {
          ++perfect_read_no_change_;
          delta_ -= new_lev_dist;
          degraded = false;
        }
        else
        {
          ++perfect_read_worsened_;
          ruined_.writeRecord(r2);
        }
      }
      else
        ++perfect_read_no_change_;
    } else {
      auto bc = r1.basesCovered();
      auto r1_seq = r1.mapped_sequence;
      auto r2_seq = r2.mapped_sequence;
      if (r2_seq.length > bc &&
          ((r2.strand == '+' && equal(r2_seq.take(bc), r2.dna().take(bc))) ||
           (r2.strand == '-' && equal(r2_seq.retro().take(bc),
                                      r2.dna().to!string.retro.take(bc)))))
      {
        delta_ -= new_lev_dist - old_lev_dist;
        degraded = false;
        ++err_read_improved_;
        ++err_read_fully_corrected_;
        if (old_lev_dist == 1)
          ++one_err_read_corrected_;
      } else if (r2_seq.length >= r1_seq.length &&
          ((r2.strand == '+' && equal(r1_seq, r2_seq.take(r1_seq.length))) ||
           (r2.strand == '-' && equal(r1_seq.retro(), 
                                      r2_seq.retro().take(r1_seq.length)))))
      {
        delta_ -= new_lev_dist - old_lev_dist;
        degraded = false;
        ++err_read_no_change_;
        if (old_lev_dist == 1)
          ++one_err_read_no_change_;
      } else if (new_lev_dist > old_lev_dist) {
        ++err_read_worsened_;
        if (old_lev_dist == 1)
          ++one_err_read_worsened_;

        damaged_.writeRecord(r2);
      } else if (new_lev_dist < old_lev_dist) {
        ++err_read_improved_;
        if (new_lev_dist == 0)
          ++err_read_fully_corrected_;

        if (old_lev_dist == 1)
          ++one_err_read_corrected_;

      } else if (equal(old_cigar_path, new_cigar_path)) {
        ++err_read_no_change_;
        if (old_lev_dist == 1)
          ++one_err_read_no_change_;

      } else {
        ++err_read_same_lev_dist_;
        if (old_lev_dist == 1)
          ++one_err_read_same_lev_dist_;
      }
    }

    if (degraded) {
      worsened_.writeln("=================================================================");
      worsened_.writeln("Read name: \n", r1.name);
      worsened_.writeln("Read sequence: \n", sequenceStr(r1));
      worsened_.writeln("CIGAR path: \n", old_cigar_path.map!(to!char)());
      worsened_.writeln("Corrected read sequence: \n", sequenceStr(r2));
      worsened_.writeln("corr. path: \n", new_cigar_path.map!(to!char)());
    }
  }

  /// Output reports
  void writeReports() {

    import std.stdio;
    writeln("Total number of reads:                   ", total_);
    writeln("-----------------------------------------");
    writeln("Reads with errors:                       ", n_err_);
    writeln("   not changed:                          ", err_read_no_change_);
    writeln("   same Lev. distance:                   ", err_read_same_lev_dist_);
    writeln("   improved:                             ", err_read_improved_);
    writeln("       fully corrected:                  ", err_read_fully_corrected_);
    writeln("       partially corrected:              ", err_read_improved_ - err_read_fully_corrected_);
    writeln("   worsened:                             ", err_read_worsened_);
    writeln("   --------------------------------------");
    writeln("   Reads with one error:                 ", n_one_err_);
    writeln("       not changed:                      ", one_err_read_no_change_);
    writeln("       same Lev. distance:               ", one_err_read_same_lev_dist_);
    writeln("       corrected:                        ", one_err_read_corrected_);
    writeln("       worsened:                         ", one_err_read_worsened_);
    writeln("-----------------------------------------");
    writeln("Perfect reads:                           ", n_perfect_);
    writeln("   not changed:                          ", perfect_read_no_change_);
    writeln("   worsened:                             ", perfect_read_worsened_);
    writeln("-----------------------------------------");
    writeln("Average change in levenshtein distance:  ", to!double(delta_) / total_,
            " / (max. possible ", to!double(old_delta_) / total_, ")");
    writeln("Average difference in # aligned bases:   ", to!double(delta_cov_) / total_);

    outputDistribution("lev_dist.tab", lev_distance_distribution_);
  }

  private {
    File worsened_;

    int[1024] lev_distance_distribution_;

    long delta_;
    long old_delta_;
    ulong total_;
    long delta_cov_;

    ulong n_err_;
    ulong err_read_no_change_;
    ulong err_read_same_lev_dist_;
    ulong err_read_improved_;
    ulong err_read_fully_corrected_;
    ulong err_read_worsened_;

    ulong n_one_err_;
    ulong one_err_read_no_change_;
    ulong one_err_read_same_lev_dist_;
    ulong one_err_read_corrected_;
    ulong one_err_read_worsened_;

    ulong n_perfect_;
    ulong perfect_read_no_change_;
    ulong perfect_read_worsened_;
  }
}

void main(string[] args) {
  if (args.length < 3) { stderr.writeln("Usage: ", args[0], " <raw.bam> <corrected.bam>"); return; }
  auto pool = new TaskPool(16);
  scope(exit) pool.finish();
  auto raw = new BamReader(args[1], pool);
  auto corrected = new BamReader(args[2], pool);

  auto comparator = new Comparator(raw, corrected);
  comparator.run();
  comparator.writeReports();
}
