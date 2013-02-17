import bio.bam.reader, bio.bam.read, bio.bam.md.reconstruct, bio.bam.baseinfo;
import bio.core.sequence, bio.core.fasta, bio.core.base;
import std.stdio, std.range, std.conv, std.algorithm, 
       std.parallelism, std.array, std.traits, std.exception;

// FIXME: more fine-grained statistics about 
//        corrected/miscorrected/uncorrected errors

struct Pair {
  BamRead read;
  string corrected_sequence;
}

struct ProcessedPair {
  BamRead read;
  string corrected_sequence;
  size_t edit_distance; // levenshtein distance b/w corrected read and reference
  EditOp[] edit_path;
}

/// Align corrected read against reference
ProcessedPair process(Pair pair) {
  auto expected = to!string(dna(pair.read));

  if (pair.read.strand == '-')
    expected = to!string(nucleotideSequence(expected).reverse);

  auto seq = pair.corrected_sequence;
  auto lev_dist_and_path = levenshteinDistanceAndPath(expected, seq);

  auto dist = lev_dist_and_path[0];
  auto path = lev_dist_and_path[1];

  while (path.back == 'r') {
      path.popBack();
      --dist;
  }

  return ProcessedPair(pair.read, pair.corrected_sequence, dist, path);
}

/// Takes BAM reads and corrected reads as input,
/// yields pairs of corresponding ones.
struct PairRange(R, CR) {
  private {
    R reads_;
    CR corrected_reads_;

    Pair front_;
    bool empty_;
  }

  this(R reads, CR corrected_reads) {
    reads_ = reads;
    corrected_reads_ = corrected_reads;
    popFront();
  }

  bool empty() @property {
    return empty_;
  }

  Pair front() @property {
    return front_;
  }

  void popFront() {
    if (corrected_reads_.empty) {
      empty_ = true;
      return;
    }

    auto record = corrected_reads_.front();
    while (!reads_.empty && reads_.front.name != record.header)
      reads_.popFront();

    if (reads_.empty) {
      empty_ = true;
      return;
    }

    front_ = Pair(reads_.front, record.sequence);
    corrected_reads_.popFront();
  }
}

class Comparator(R, CR) 
  if (is(ElementType!R == BamRead) && is(ElementType!CR == FastaRecord))
{
  private {
    PairRange!(R, CR) read_pairs_;
  }

  this(R reads, CR corrected_reads) {
    read_pairs_ = PairRange!(R, CR)(reads, corrected_reads);
    worsened_ = File("worsened.txt", "w+");
  }

  /// Collect statistics
  void run() {
    version (serial) {
      auto processed_pairs = map!process(read_pairs_);
    } else {
      auto taskpool = new TaskPool(totalCPUs);
      scope(exit) taskpool.finish();

      auto processed_pairs = taskpool.map!process(read_pairs_);
    }

    foreach (pair; processed_pairs) {
      updateStats(pair);
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

  private static size_t pathToDistance(in EditOp[] path) {
    size_t dist;
    foreach (edit_op; path)
      if (edit_op != 'n')
        ++dist;
    return dist;
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

  // same for corrected sequence with known edit path
  private static string corrSequenceStr(string corr, in EditOp[] path) {
    string res;
    foreach (op; path) {
      if (op != 'r') {
        res ~= corr.front;
        corr.popFront();
      } else {
        res ~= '*';
      }
    }
    return res;
  }

  private void updateStats(ProcessedPair pair) {
    auto lev_dist = pair.edit_distance;
    auto path = pair.edit_path;

    lev_distance_distribution_[lev_dist] += 1;

    // ------------- compare with CIGAR ----------------------------------
    auto cigar_path = cigarPath(pair.read);
    auto old_lev_dist = pathToDistance(cigar_path);

    ++total_;
    delta_ += to!long(lev_dist) - to!long(old_lev_dist);
    old_delta_ -= old_lev_dist;

    if (old_lev_dist == 0)
      ++n_perfect_;
    else
      ++n_err_;

    if (old_lev_dist == 1)
      ++n_one_err_;

    if (old_lev_dist == 0) {
      if (lev_dist > 0)
        ++perfect_read_worsened_;
      else
        ++perfect_read_no_change_;
    } else {
      if (lev_dist > old_lev_dist) {
        ++err_read_worsened_;
        if (old_lev_dist == 1)
          ++one_err_read_worsened_;

      } else if (lev_dist < old_lev_dist) {
        ++err_read_improved_;
        if (lev_dist == 0)
          ++err_read_fully_corrected_;

        if (old_lev_dist == 1)
          ++one_err_read_corrected_;

      } else if (equal(cigar_path, path)) {
        ++err_read_no_change_;
        if (old_lev_dist == 1)
          ++one_err_read_no_change_;

      } else {
        ++err_read_same_lev_dist_;
        if (old_lev_dist == 1)
          ++one_err_read_same_lev_dist_;
      }
    }

    if (old_lev_dist < lev_dist) {
      worsened_.writeln("=================================================================");
      worsened_.writeln("Read name: \n", pair.read.name);
      worsened_.writeln("Read sequence: \n", sequenceStr(pair.read));
      worsened_.writeln("CIGAR path: \n", cigar_path.map!(to!char)());
      worsened_.writeln("Corrected read sequence: \n", 
                        corrSequenceStr(pair.corrected_sequence, path));
      worsened_.writeln("corr. path: \n", path.map!(to!char)());
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

    outputDistribution("lev_dist.tab", lev_distance_distribution_);
  }

  private {
    File worsened_;

    int[1024] lev_distance_distribution_;

    long delta_;
    long old_delta_;
    ulong total_;

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

auto makeComparator(R, CR)(R reads, CR corrected_reads) {
  return new Comparator!(R, CR)(reads, corrected_reads);
}

void main(string[] args) {
  auto bam = new BamReader(args[1]);
  auto corrected = fastaRecords(args[2]);

  auto comparator = makeComparator(bam.reads, corrected);
  comparator.run();
  comparator.writeReports();
}
