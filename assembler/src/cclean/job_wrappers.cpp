#include <set>
#include "job_wrappers.hpp"
#include "logger/log_writers.hpp"

bool is_alignment_good(const StripedSmithWaterman::Alignment& a, const std::string& query, double aligned_part_fraction) {
  // FIXME: Check for the end of reference here (somehow)
  return (std::min(a.query_end - a.query_begin + 1, a.ref_end - a.ref_begin + 1) / (double) query.size() > aligned_part_fraction);
}

bool AlignmentJobWrapper::operator()(const Read &r) {
  try {
    const std::string &name = r.getName();
    const std::string &ref = r.getSequenceString();
    auto it = data->get_data_iterator();
    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;
    aligner.SetReferenceSequence(ref.c_str(), ref.size());
    for (unsigned i = 0; i < data->get_sequences_amount(); ++i) {
      StripedSmithWaterman::Alignment alignment;
      const std::string &query = *(it->second);
      aligner.Align(query.c_str(), filter, &alignment);
      std::string& database_name = *(it->first);

      if (alignment.mismatches < mismatch_threshold && is_alignment_good(alignment, query, aligned_part_fraction)) {
#       pragma omp critical
        {
          aligned_ += 1;
          print_alignment(output, alignment, ref, query, name, database_name);
          print_bed(bed, name, alignment.ref_begin, alignment.ref_end);
        }
      }
      it++;
    }
  } catch (std::exception& e) {
    ERROR(e.what());
  }
  return false;
}

bool ExactMatchJobWrapper::operator()(const Read &r) {
  try {
    const std::string& name = r.getName();
    const std::string& sequence = r.getSequenceString();
    seq2index_t res = ahoCorasick.search(sequence);

    if (res.size() > 0) {
#     pragma omp critical
      {
        aligned_ += 1;
        print_match(output, bed, res, name, sequence, data);
      }
    }
  } catch (std::exception& e) {
    ERROR(e.what());
  }
  return false;
}

bool ExactAndAlignJobWrapper::operator()(const Read &r) {
  try {
    const std::string& name = r.getName();
    const std::string& sequence = r.getSequenceString();

    //try to exact match the sequences
    auto matchingSequences = dbAhoCorasick.search(sequence);
    if (!matchingSequences.empty()) {
#     pragma omp critical
      {
        aligned_ += 1;
        print_match(output, bed, matchingSequences, name, sequence, data);
      }
      return false; //exact match is better than an alignment -> no need to align
    }

    //try to search in kmers db
    auto matchingKmers = kmersAhoCorasick.search(sequence);
    std::set<std::string * , Compare> setOfContaminations2check;
    for (auto it = matchingKmers.begin(); it != matchingKmers.end(); ++it) {
      std::set<std::string *, Compare> setOfSeqs;
      data->get_sequences_for_kmer(*(it->first), setOfSeqs);
      setOfContaminations2check.insert(setOfSeqs.begin(), setOfSeqs.end());
    }

    // Try to align the contaminations for corresponding kmers
    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;
    aligner.SetReferenceSequence(sequence.c_str(), sequence.size());
    for (auto it = setOfContaminations2check.begin(); it != setOfContaminations2check.end(); ++it) {
      StripedSmithWaterman::Alignment alignment;
      const std::string& query = *(*it);
      aligner.Align(query.c_str(), filter, &alignment);

      std::string database_name;
      data->get_name_by_sequence(query, database_name);

      if (alignment.mismatches < mismatch_threshold && is_alignment_good(alignment, query, aligned_part_fraction)) {
#       pragma omp critical
        {
          aligned_ += 1;
          print_alignment(output, alignment, sequence, query, name, database_name);
          print_bed(bed, name, alignment.ref_begin, alignment.ref_end);
        }
      }
    }

  } catch (std::exception& e) {
    ERROR(e.what());
  }
  return false;
}
