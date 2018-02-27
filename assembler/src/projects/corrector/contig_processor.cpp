//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "contig_processor.hpp"
#include "config_struct.hpp"
#include "variants_table.hpp"

#include "io/reads/ireader.hpp"
#include "io/reads/osequencestream.hpp"
#include "io/reads/file_reader.hpp"
#include "io/reads/single_read.hpp"
#include "utils/filesystem/path_helper.hpp"

#include <boost/algorithm/string.hpp>

using namespace std;

namespace corrector {

void ContigProcessor::ReadContig() {
    io::FileReadStream frs(contig_file_);
    io::SingleRead cur_read;
    frs >> cur_read;
    if (!frs.eof()) {
#pragma omp critical
        {
            ERROR("Non unique sequnce in one contig fasta!");
        }
    }
    contig_name_ = cur_read.name();
    contig_ = cur_read.GetSequenceString();

    output_contig_file_ = fs::append_path(fs::parent_path(contig_file_), fs::basename(contig_file_) + ".ref.fasta");
    charts_.resize(contig_.length());
}

void ContigProcessor::UpdateOneRead(const SingleSamRead &tmp, MappedSamStream &sm) {
    unordered_map<size_t, position_description> all_positions;
    if (tmp.contig_id() < 0) {
        return;
    }
    auto cur_s = sm.get_contig_name(tmp.contig_id());
    if (contig_name_.compare(cur_s) != 0) {
        return;
    }
    CountPositions(tmp, all_positions);
    size_t error_num = 0;

    for (auto &pos : all_positions) {
        charts_[pos.first].update(pos.second);
        if (pos.second.FoundOptimal(contig_[pos.first]) != var_to_pos[(int) contig_[pos.first]]) {
            error_num++;
        }
    }

    if (error_num >= error_counts_.size())
        error_counts_[error_counts_.size() - 1]++;
    else
        error_counts_[error_num]++;
}

//returns: number of changed nucleotides;
size_t ContigProcessor::UpdateOneBase(size_t i, stringstream &ss, const unordered_map<size_t, position_description> &interesting_positions) const{
    char old = (char) toupper(contig_[i]);
    auto strat = corr_cfg::get().strat;
    size_t maxi = charts_[i].FoundOptimal(contig_[i]);
    auto i_position = interesting_positions.find(i);
    if (i_position != interesting_positions.end()) {
        size_t maxj = i_position->second.FoundOptimal(contig_[i]);
        if (maxj != maxi) {
            DEBUG("Interesting positions differ with majority!");
            DEBUG("On position " << i << "  old: " << old << " majority: " << pos_to_var[maxi] << "interesting: " << pos_to_var[maxj]);
            if (strat != Strategy::MajorityOnly) {
                if (charts_[i].votes[maxj] > interesting_weight_cutoff)
                    maxi = maxj;
                else
                    DEBUG(" alternative interesting position with weight " << charts_[i].votes[maxj] <<
                          " fails weight cutoff");
            }
        }
    }
    if (old != pos_to_var[maxi]) {
        DEBUG("On position " << i << " changing " << old << " to " << pos_to_var[maxi]);
        DEBUG(charts_[i].str());
        if (maxi < Variants::Deletion) {
            ss << pos_to_var[maxi];
            return 1;
        } else if (maxi == Variants::Deletion) {
            return 1;
        } else if (maxi == Variants::Insertion) {
            string maxj = "";
            //first base before insertion;
            size_t new_maxi = var_to_pos[(int) contig_[i]];
            int new_maxx = charts_[i].votes[new_maxi];
            for (size_t k = 0; k < MAX_VARIANTS; k++) {
                if (new_maxx < charts_[i].votes[k] && (k != Variants::Insertion) && (k != Variants::Deletion)) {
                    new_maxx = charts_[i].votes[k];
                    new_maxi = k;
                }
            }
            ss << pos_to_var[new_maxi];
            int max_ins = 0;
            for (const auto &ic : charts_[i].insertions) {
                if (ic.second > max_ins) {
                    max_ins = ic.second;
                    maxj = ic.first;
                }
            }
            DEBUG("most popular insertion: " << maxj);
            ss << maxj;
            if (old == maxj[0]) {
                return (int) maxj.length() - 1;
            } else {
                return (int) maxj.length();
            }
        } else {
            //something strange happened
            WARN("While processing base " << i << " unknown decision was made");
            return 0;
        }
    } else {
        ss << old;
        return 0;
    }
}


bool ContigProcessor::CountPositions(const SingleSamRead &read, unordered_map<size_t, position_description> &ps) const {

    if (read.contig_id() < 0) {
        DEBUG("not this contig");
        return false;
    }
    //TODO: maybe change to read.is_properly_aligned() ?
    if (read.map_qual() == 0) {
        DEBUG("zero qual");
        return false;
    }
    int pos = read.pos();
    if (pos < 0) {
        WARN("Negative position " << pos << " found on read " << read.name() << ", skipping");
        return false;
    }
    size_t position = size_t(pos);
    int mate = 1;  // bonus for mate mapped can be here;
    size_t l_read = (size_t) read.data_len();
    size_t l_cigar = read.cigar_len();

    int aligned_length = 0;
    uint32_t *cigar = read.cigar_ptr();
    //* in cigar;
    if (l_cigar == 0)
        return false;
    if (bam_cigar_opchr(cigar[0]) == '*')
        return false;
    for (size_t i = 0; i < l_cigar; i++)
        if (bam_cigar_opchr(cigar[i]) == 'M')
            aligned_length += bam_cigar_oplen(cigar[i]);
//It's about bad aligned reads, but whether it is necessary?
    double read_len_double = (double) l_read;
    if ((aligned_length < min(read_len_double * 0.4, 40.0)) && (position > read_len_double / 2) && (contig_.length() > read_len_double / 2 + (double) position)) {
        return false;
    }
    int state_pos = 0;
    int shift = 0;
    size_t skipped = 0;
    size_t deleted = 0;
    string insertion_string = "";
    auto seq = read.seq_ptr();
    for (size_t i = 0; i < l_read; i++) {
        DEBUG(i << " " << position << " " << skipped);
        if (shift + bam_cigar_oplen(cigar[state_pos]) <= i) {
            shift += bam_cigar_oplen(cigar[state_pos]);
            state_pos += 1;
        }
        if (insertion_string != "" and bam_cigar_opchr(cigar[state_pos]) != 'I') {
            VERIFY(i + position >= skipped + 1);
            size_t ind = i + position - skipped - 1;
            if (ind >= contig_.length())
                break;
            ps[ind].insertions[insertion_string] += 1;
            insertion_string = "";
        }
        char cur_state = bam_cigar_opchr(cigar[state_pos]);
        if (cur_state == 'M') {
            VERIFY(i >= deleted);
            if (i + position < skipped) {
                WARN(i << " " << position << " " << skipped);
                INFO(read.name());
            }
            VERIFY(i + position >= skipped);

            size_t ind = i + position - skipped;
            size_t cur = var_to_pos[(int) bam_nt16_rev_table[bam1_seqi(seq, i - deleted)]];
            if (ind >= contig_.length())
                continue;
            ps[ind].votes[cur] = ps[ind].votes[cur] + mate;

        } else {
            if (cur_state == 'I' || cur_state == 'H' || cur_state == 'S' ) {
                if (cur_state == 'I') {
                    if (insertion_string == "") {
                        size_t ind = i + position - skipped - 1;
                        if (ind >= contig_.length())
                            break;
                        ps[ind].votes[Variants::Insertion] += mate;
                    }
                    insertion_string += bam_nt16_rev_table[bam1_seqi(seq, i - deleted)];
                }
                skipped += 1;
            } else if (bam_cigar_opchr(cigar[state_pos]) == 'D') {
                if (i + position - skipped >= contig_.length())
                    break;
                ps[i + position - skipped].votes[Variants::Deletion] += mate;
                deleted += 1;
            }
        }
    }
    if (insertion_string != "" and bam_cigar_opchr(cigar[state_pos]) != 'I') {
        VERIFY(l_read + position >= skipped + 1);
        size_t ind = l_read + position - skipped - 1;
        if (ind < contig_.length()) {
            ps[ind].insertions[insertion_string] += 1;
        }
        insertion_string = "";
    }
    return true;
}


bool ContigProcessor::CountPositions(const PairedSamRead &read, unordered_map<size_t, position_description> &ps) const {

    TRACE("starting pairing");
    bool t1 = CountPositions(read.Left(), ps );
    unordered_map<size_t, position_description> tmp;
    bool t2 = CountPositions(read.Right(), tmp);
    //overlaps.. multimap? Look on qual?
    if (ps.size() == 0 || tmp.size() == 0) {
        //We do not need paired reads which are not really paired
        ps.clear();
        return false;
    }
    TRACE("counted, uniting maps of " << tmp.size() << " and " << ps.size());
    ps.insert(tmp.begin(), tmp.end());
    TRACE("united");
    return (t1 && t2);
}

size_t ContigProcessor::ProcessMultipleSamFiles() {
    error_counts_.resize(kMaxErrorNum);
    for (const auto &sf : sam_files_) {
        MappedSamStream sm(sf.first);
        while (!sm.eof()) {
            SingleSamRead tmp;
            sm >> tmp;

            UpdateOneRead(tmp, sm);
        }
        sm.close();
    }
    size_t total_coverage = 0;
    for (const auto &pos: charts_)
        total_coverage += pos.TotalMapped();
    size_t average_coverage = total_coverage / contig_.length();
    size_t different_cov = 0;
    for (const auto &pos: charts_)
        if ((pos.TotalMapped() < average_coverage / 2) || (pos.TotalMapped() > (average_coverage * 3) / 2))
            different_cov++;
    if (different_cov < contig_.length() * 3/ 10) {
        interesting_weight_cutoff = int (average_coverage / 2);
        DEBUG ("coverage is relatively uniform, average coverage is " << average_coverage
               << " setting interesting positions heuristics to " << interesting_weight_cutoff);
    }
    ipp_.FillInterestingPositions(charts_);
    for (const auto &sf : sam_files_) {
        MappedSamStream sm(sf.first);
        while (!sm.eof()) {
            unordered_map<size_t, position_description> ps;
            if (sf.second == io::LibraryType::PairedEnd ) {
                PairedSamRead tmp;
                sm >> tmp;
                CountPositions(tmp, ps);
            } else {
                SingleSamRead tmp;
                sm >> tmp;
                CountPositions(tmp, ps);
            }
            ipp_.UpdateInterestingRead(ps);
        }
        sm.close();
    }
    ipp_.UpdateInterestingPositions();
    unordered_map<size_t, position_description> interesting_positions = ipp_.get_weights();
    stringstream s_new_contig;
    size_t total_changes = 0;
    for (size_t i = 0; i < contig_.length(); i++) {
        total_changes += UpdateOneBase(i, s_new_contig, interesting_positions);
    }
    vector<string> contig_name_splitted;
    boost::split(contig_name_splitted, contig_name_, boost::is_any_of("_"));
    io::OFastaReadStream oss(output_contig_file_);
    for(size_t i = 0; i < contig_name_splitted.size(); i++) {
        if (contig_name_splitted[i] == "length" && i + 1 < contig_name_splitted.size()) {
            contig_name_splitted[i + 1] = std::to_string(int(s_new_contig.str().length()));
            break;
        }
    }
    std::string new_header = contig_name_splitted[0];
    for(size_t i = 1; i < contig_name_splitted.size(); i++) {
        new_header += "_" + contig_name_splitted[i];
    }
    oss << io::SingleRead(new_header, s_new_contig.str());

    return total_changes;
}

}
;
