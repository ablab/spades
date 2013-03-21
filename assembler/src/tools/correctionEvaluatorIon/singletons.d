// Usage:
// 1. Uncomment code for FASTA output of centers in hammer-it/main.cpp
// 2. Using TMAP, align the file with centers against the reference.
//    Don't forget to specify BAM output. ('-o 1')
//    It is also recommended to switch off soft-clipping with '-g 3'
// 3. Run ./singletons on the BAM file with aligned centers.

import bio.bam.reader;
import std.parallelism, std.stdio, std.conv, std.algorithm, std.range;

void main(string[] args) {
  auto pool = new TaskPool(4);
  scope(exit) pool.finish();

  uint centers;
  uint genomic_centers;
  uint unmapped_centers;
  uint clipped_centers;
  uint[32] center_distances;

  uint singletons;
  uint genomic_singletons;
  uint unmapped_singletons;
  uint clipped_singletons;
  uint[32] singleton_distances;

  auto bam = new BamReader(args[1], pool);
  foreach (read; bam.reads) {
    auto distance = read["NM"];
    auto has_distance = !distance.is_nothing;
    auto is_genomic = !read.is_unmapped && distance == 0 && read.cigar.length == 1;
    auto is_singleton = read.name[$ - 1] == 'n';
    auto is_clipped = !read.cigar.filter!(op => op.type == 'S').empty();
    centers++;
    if (is_genomic)
      genomic_centers++;
    else if (is_clipped)
      clipped_centers++;
    else if (read.is_unmapped)
      unmapped_centers++;
    else if (has_distance)
      center_distances[distance.to!int()]++;

    if (is_singleton) {
      singletons++;
      if (is_genomic)
        genomic_singletons++;
      else if (is_clipped)
        clipped_singletons++;
      else if (read.is_unmapped)
        unmapped_singletons++;
      else if (has_distance)
        singleton_distances[distance.to!int()]++;
    }
  }

  writeln("centers:            ", centers);

  /*
  writeln("  genomic:          ", genomic_centers);
  writeln("  non-genomic:      ", centers - genomic_centers);
  writeln("    unmapped:       ", unmapped_centers);
  writeln("    mapped:         ", centers - unmapped_centers);
  if (clipped_centers > 0)
    writeln("      clipped:      ", clipped_centers);

  foreach (size_t dist, count; center_distances[])
    if (count > 0)
      writeln("      distance ", dist, ": ", count);
  */

  writeln();
  writeln("  singletons:       ", singletons);
  writeln("    genomic:        ", genomic_singletons);
  writeln("    non-genomic:    ", singletons - genomic_singletons);
  writeln("      unmapped:     ", unmapped_singletons);
  writeln("      mapped:       ", singletons - genomic_singletons - unmapped_singletons);
  if (clipped_singletons > 0)
    writeln("        clipped:    ", clipped_singletons);

  foreach (size_t dist, count; singleton_distances[])
    if (count > 0)
      writeln("        distance ", dist, ": ", count);

  writeln();
  writeln(" !singletons:       ", centers - singletons);
  writeln("    genomic:        ", genomic_centers - genomic_singletons);
  writeln("    non-genomic:    ", (centers - genomic_centers) - (singletons - genomic_singletons));
  writeln("      unmapped:     ", unmapped_centers - unmapped_singletons);
  writeln("      mapped:       ", (centers - genomic_centers - unmapped_centers) - (singletons - genomic_singletons - unmapped_singletons));
  if (clipped_centers > 0)
    writeln("        clipped:    ", clipped_centers - clipped_singletons);

  foreach (size_t dist, count; center_distances[])
    if (count > singleton_distances[dist])
      writeln("        distance ", dist, ": ", count - singleton_distances[dist]);
}
