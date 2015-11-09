//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "bwa_pair_info_filler.hpp"


namespace bwa_pair_info {


void MapperReadT::ParseCigar(const string& cigar) {
    string num = "";
    bool left_side = true;
    for (size_t i = 0; i < cigar.length(); ++i) {
        if (isdigit(cigar[i])) {
            num += cigar[i];
        }
        else {
            if (cigar[i] == 'H') {
                if (left_side)
                    left_hard_clip_ = (uint16_t) std::stoi(num);
                else
                    right_hard_clip_ = (uint16_t) std::stoi(num);
                num = "";
            }
            else if (cigar[i] == 'S') {
                if (left_side)
                    left_soft_clip_ = (uint16_t) std::stoi(num);
                else
                    right_soft_clip_ = (uint16_t) std::stoi(num);
                num = "";
            }
            else {
                left_side = false;
                num = "";
            }
        }
    }
}

void BWACorrectingProcessor::ProcessPairedRead(const MapperReadT& l, const MapperReadT& r) {
    using io::LibraryOrientation;

    if (!l.IsValid() || !r.IsValid()) {
        return;
    }
    ++count_;

    MappedPositionT left_pos(edge_id_map_.at(stoi(l.get_contig_id())), l.get_pos());
    MappedPositionT right_pos(edge_id_map_.at(stoi(r.get_contig_id())), r.get_pos());

    if (!CheckAlignments(left_pos, right_pos)) {
        return;
    }

    int r_from_pos_to_right_end = r.get_len() + r.get_right_hard_clip() - r.get_left_soft_clip();
    int l_from_pos_to_left_end = l.get_left_soft_clip() + l.get_left_hard_clip();

    if ((!l.is_forward() && (lib_.orientation() == LibraryOrientation::FF || lib_.orientation() == LibraryOrientation::FR)) ||
        (l.is_forward() && (lib_.orientation() == LibraryOrientation::RF || lib_.orientation() == LibraryOrientation::RR))) {
        left_pos.e_ = g_.conjugate(left_pos.e_);
        left_pos.pos_ = (int) g_.length(left_pos.e_) - left_pos.pos_ - (l.get_len() - l.get_left_soft_clip() - l.get_right_soft_clip()) + (int) g_.k();
        l_from_pos_to_left_end = l.get_right_soft_clip() + l.get_right_hard_clip();
    }
    if ((!r.is_forward() && (lib_.orientation() == LibraryOrientation::FF || lib_.orientation() == LibraryOrientation::RF)) ||
        (r.is_forward() && (lib_.orientation() == LibraryOrientation::FR || lib_.orientation() == LibraryOrientation::RR))) {
        right_pos.e_ = g_.conjugate(right_pos.e_);
        right_pos.pos_ = (int) g_.length(right_pos.e_) - right_pos.pos_ - (r.get_len() - r.get_left_soft_clip() - r.get_right_soft_clip()) + (int) g_.k();
        r_from_pos_to_right_end = r.get_len() + r.get_left_hard_clip() - r.get_right_soft_clip();
    }

    right_pos.pos_ = right_pos.pos_ + r_from_pos_to_right_end;
    left_pos.pos_ = left_pos.pos_ - l_from_pos_to_left_end;

    ProcessAlignments(left_pos, right_pos);
}

EdgePair BWAIndexFiller::ConjugatePair(EdgePair ep) const {
    return make_pair(g_.conjugate(ep.second), g_.conjugate(ep.first));
}

void BWAIndexFiller::ProcessAlignments(const MappedPositionT& l, const MappedPositionT& r) {
    EdgePair ep{l.e_, r.e_};
    TRACE("Lpos " << l.pos_ << ", Rpos " << r.pos_);
    int edge_distance = (int) lib_.data().mean_insert_size  - r.pos_ + l.pos_;
    TRACE("Distance " << edge_distance);

    if (ep > ConjugatePair(ep)) {
        ep = ConjugatePair(ep);
        edge_distance = edge_distance + (int) g_.length(r.e_) - (int) g_.length(l.e_);
        TRACE("New distance " << edge_distance);
    }

    paired_index_.AddPairInfo(ep.first, ep.second, { (double) edge_distance, 1.0 });
}

bool BWAIndexFiller::CheckAlignments(const MappedPositionT& l, const MappedPositionT& r) {
    return g_.length(l.e_) >= min_contig_len_ && g_.length(r.e_) >= min_contig_len_;
}


void BWAPairInfoFiller::OutputEdges(const string &filename) const {
    io::osequencestream_simple oss(filename);
    for (auto it = g_.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
        debruijn_graph::EdgeId e = *it;
        oss.set_header(ToString(g_.int_id(e)));
        oss << g_.EdgeNucls(e);
    }
}
void BWAPairInfoFiller::FillEdgeIdMap() {
    for (auto it = g_.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
        debruijn_graph::EdgeId e = *it;
        edge_id_map_.insert(make_pair(g_.int_id(e), e));
    }
}

bool BWAPairInfoFiller::CreateIndex(const string& contigs) {
    int run_res = 0;
    string err_log = path::append_path(work_dir_, "index.err");
    string index_line = bwa_path_ + string(" index ") + "-a is " + contigs + " 2>" + err_log;
    INFO("Running bwa index ... ");
    INFO("Command line: " << index_line);
    run_res = system(index_line.c_str());
    if (run_res != 0) {
        ERROR("bwa index failed, cannot align reads");
        return false;
    }
    return true;
}


bool BWAPairInfoFiller::RunBWA(const string& reads_file, const string& out_sam_file) const {
    string run_command = bwa_path_ + " mem -t " + ToString(nthreads_) + " " + index_base_ + " "  + reads_file + "  > " + out_sam_file + " 2>"
        + out_sam_file + ".txt";
    INFO("Running bwa mem ...");
    INFO("Command line: " << run_command);

    int run_res = system(run_command.c_str());
    if (run_res != 0) {
        ERROR("bwa mem failed, cannot align reads");
        return false;
    }
    return true;
}

bool BWAPairInfoFiller::AlignLib(const SequencingLibraryT& lib,
                                 const string& sam_file_base,
                                 vector<pair<string, string>>& resulting_sam_files) {

    VERIFY_MSG(Init(), "BWA index was not constructed properly");
    resulting_sam_files.clear();
    size_t file_index = 0;
    bool any_aligned = false;

    for (auto iter = lib.paired_begin(); iter != lib.paired_end(); iter++) {
        string left_reads = iter->first;
        string left_sam = sam_file_base + "_1_" + ToString(file_index) + ".sam";
        bool res = RunBWA(left_reads, left_sam);
        if (!res) {
            WARN("Failed to align left reads " << left_reads);
            continue;
        }
        string right_reads = iter->second;
        string right_sam = sam_file_base + "_2_" + ToString(file_index) + ".sam";
        res = RunBWA(right_reads, right_sam);
        if (!res) {
            WARN("Failed to align right reads " << right_reads);
            continue;
        }

        resulting_sam_files.push_back(make_pair(left_sam, right_sam));
        any_aligned = true;
    }
    return any_aligned;
}


void BWAPairInfoFiller::ProcessSAMFiles(const string &left_sam, const string &right_sam,
                                        BWAPairedReadProcessor& processor) {
    unordered_map<string, MapperReadT> left_reads;
    unordered_map<string, MapperReadT> right_reads;

    INFO("Reading SAM files " << left_sam << " and " << right_sam);
    MappedSamStream lf(left_sam);
    MappedSamStream rf(right_sam);
    while (!lf.eof() || !rf.eof()) {
        SingleSamRead left_read;
        MapperReadT left_data;
        string l_name = "";

        SingleSamRead right_read;
        MapperReadT right_data;
        string r_name = "";

        if (!lf.eof()) {
            lf >> left_read;
            l_name = left_read.get_name();
            if (left_read.is_properly_aligned()) {
                TRACE("Left read " << l_name);
                left_data = MapperReadT(string(lf.get_contig_name(left_read.get_contig_id())),
                                        left_read.get_pos(),
                                        left_read.get_data_len(),
                                        left_read.get_strand(),
                                        left_read.get_cigar());
            }
            else if (!left_read.is_main_alignment()) {
                TRACE("Ignoring left read");
                l_name = "";
            }
        }
        if (!rf.eof()) {
            rf >> right_read;
            r_name = right_read.get_name();
            if (right_read.is_properly_aligned()) {
                TRACE("Right read " << r_name);
                right_data = MapperReadT(string(rf.get_contig_name(right_read.get_contig_id())),
                                         right_read.get_pos(),
                                         right_read.get_data_len(),
                                         right_read.get_strand(),
                                         right_read.get_cigar());
            }
            else if (!right_read.is_main_alignment()) {
                TRACE("Ignoring right read");
                r_name = "";
            }
        }

        if (l_name == r_name) {
            TRACE("Equal processing");
            processor.ProcessPairedRead(left_data, right_data);
            continue;
        }

        if (r_name != "") {
            auto it = left_reads.find(r_name);
            if (it != left_reads.end())  {
                TRACE("Right read's mate found, processing");
                processor.ProcessPairedRead(it->second, right_data);
                left_reads.erase(it);
            }
            else {
                TRACE("Right read's mate not found, adding to map");
                if (right_reads.count(r_name) == 0) {
                    right_reads.emplace(r_name, right_data);
                } else {
                    WARN("Right read " << r_name << " is duplicated!");
                }
            }
        }

        if (l_name != "") {
            auto it = right_reads.find(l_name);
            if (it != right_reads.end()) {
                TRACE("Left read's mate found, processing");
                processor.ProcessPairedRead(left_data, it->second);
                right_reads.erase(it);
            }
            else {
                TRACE("Left read's mate not found, adding to map");
                if (left_reads.count(l_name) == 0) {
                    left_reads.emplace(l_name, left_data);
                } else {
                    WARN("Left read " << r_name << " is duplicated!");
                }

            }
        }
    }
}

bool BWAPairInfoFiller::Init() {
    if (!index_constructed_) {
        INFO("Initializing bwa pair info counter, working dir " << work_dir_);
        path::make_dir(work_dir_);
        index_base_= path::append_path(work_dir_, "long_edges.fasta");
        INFO("Saving edges to " << index_base_);
        OutputEdges(index_base_);
        FillEdgeIdMap();
        index_constructed_ = CreateIndex(index_base_);
    }
    return index_constructed_;
}

bool BWAPairInfoFiller::ProcessLib(size_t lib_index,
                                   SequencingLibraryT& lib,
                                   PairedInfoIndexT& paired_index,
                                   size_t counter_edge_len,
                                   size_t index_filler_edge_len) {

    Init();
    string lib_dir =  path::append_path(work_dir_, ToString(lib_index));
    path::make_dir(lib_dir);
    vector<pair<string, string>> sam_files;
    bool result = false;

    INFO("Processing lib #" << lib_index);
    if (!AlignLib(lib, path::append_path(lib_dir, "single"), sam_files)) {
        WARN("Failed to align lib #" << lib_index);
        return false;
    }

    INFO("Estimating insert size for library #" << lib_index);
    BWAISCounter counter(lib, edge_id_map_, g_, counter_edge_len);
    for (auto sam_pair : sam_files) {
        ProcessSAMFiles(sam_pair.first, sam_pair.second, counter);
    }

    if (!counter.RefineInsertSize(lib)) {
        lib.data().mean_insert_size = 0.0;
        WARN("Unable to estimate insert size paired library #" << lib_index);
    }
    else {
        INFO("  Estimated insert size for paired library #" << lib_index);
        INFO("  Insert size = " << lib.data().mean_insert_size <<
            ", deviation = " << lib.data().insert_size_deviation <<
            ", left quantile = " << lib.data().insert_size_left_quantile <<
            ", right quantile = " << lib.data().insert_size_right_quantile <<
            ", read length = " << lib.data().read_length);

        INFO("Collecting paired information for library #" << lib_index);
        for (auto it = g_.ConstEdgeBegin(); !it.IsEnd(); ++it)
            paired_index.AddPairInfo(*it, *it, {0., 0.});

        BWAIndexFiller filler(lib, edge_id_map_, g_, paired_index, index_filler_edge_len);
        for (auto sam_pair : sam_files) {
            ProcessSAMFiles(sam_pair.first, sam_pair.second, filler);
        }
        result = true;
    }

    path::remove_dir(lib_dir);
    return result;
}

bool BWAISCounter::CheckAlignments(const MappedPositionT& l, const MappedPositionT& r) {
    return l.e_ == r.e_ && g_.length(l.e_) >= min_contig_len_;
}

void BWAISCounter::ProcessAlignments(const MappedPositionT& l, const MappedPositionT& r) {
    ++mapped_count_;

    int is = r.pos_ - l.pos_;
    if (is > 0 || !ignore_negative_) {
        hist_[is] += 1;
    } else {
        ++negative_count_;
    }
}

bool BWAISCounter::RefineInsertSize(SequencingLibraryT& reads) const {
    using namespace omnigraph;
    size_t correctly_mapped = mapped_count_ - negative_count_;
    INFO(correctly_mapped << " paired reads (" << ((double) correctly_mapped * 100.0 / (double) count_) << "% of all) aligned to long edges");

    if (negative_count_ > 3 * correctly_mapped)
        WARN("Too much reads aligned with negative insert size. Is the library orientation set properly?");
    if (mapped_count_ == 0)
        return false;

    std::map<size_t, size_t> percentiles;
    find_mean(hist_, reads.data().mean_insert_size, reads.data().insert_size_deviation, percentiles);
    find_median(hist_, reads.data().median_insert_size, reads.data().insert_size_mad, reads.data().insert_size_distribution);
    if (reads.data().median_insert_size < reads.data().read_length) {
        return false;
    }

    std::tie(reads.data().insert_size_left_quantile, reads.data().insert_size_right_quantile) =
        GetISInterval(0.8, reads.data().insert_size_distribution);

    return !reads.data().insert_size_distribution.empty();
}
}
