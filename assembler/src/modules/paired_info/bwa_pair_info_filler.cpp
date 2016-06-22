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

//Correct read algnment according to orientation and clippings
void BWACorrectingProcessor::ProcessPairedRead(const MapperReadT& l, const MapperReadT& r) {
    using io::LibraryOrientation;

    if (!l.IsValid() || !r.IsValid()) {
        return;
    }
    ++count_;

    MappedPositionT left_pos(edge_id_map_.at(stoi(l.get_contig_id())), l.pos());
    MappedPositionT right_pos(edge_id_map_.at(stoi(r.get_contig_id())), r.pos());

    //This function if overloaded in BWAISCounter and BWAIndexFiller
    if (!CheckAlignments(left_pos, right_pos)) {
        return;
    }

    int r_from_pos_to_right_end = r.len() + r.right_hard_clip() - r.left_soft_clip();
    int l_from_pos_to_left_end = l.left_soft_clip() + l.left_hard_clip();

    if ((!l.is_forward() && (lib_.orientation() == LibraryOrientation::FF || lib_.orientation() == LibraryOrientation::FR)) ||
        (l.is_forward() && (lib_.orientation() == LibraryOrientation::RF || lib_.orientation() == LibraryOrientation::RR))) {
        left_pos.e = g_.conjugate(left_pos.e);
        left_pos.pos = (int) g_.length(left_pos.e) - left_pos.pos - (l.len() - l.left_soft_clip() - l.right_soft_clip()) + (int) g_.k();
        l_from_pos_to_left_end = l.right_soft_clip() + l.right_hard_clip();
    }
    if ((!r.is_forward() && (lib_.orientation() == LibraryOrientation::FF || lib_.orientation() == LibraryOrientation::RF)) ||
        (r.is_forward() && (lib_.orientation() == LibraryOrientation::FR || lib_.orientation() == LibraryOrientation::RR))) {
        right_pos.e = g_.conjugate(right_pos.e);
        right_pos.pos = (int) g_.length(right_pos.e) - right_pos.pos - (r.len() - r.left_soft_clip() - r.right_soft_clip()) + (int) g_.k();
        r_from_pos_to_right_end = r.len() + r.left_hard_clip() - r.right_soft_clip();
    }

    right_pos.pos = right_pos.pos + r_from_pos_to_right_end;
    left_pos.pos = left_pos.pos - l_from_pos_to_left_end;

    //This function if overloaded in BWAISCounter and BWAIndexFiller
    ProcessAlignments(left_pos, right_pos);
}

// ==== insert size counter overloads ====
bool BWAISCounter::CheckAlignments(const MappedPositionT& l, const MappedPositionT& r) {
    return l.e == r.e && g_.length(l.e) >= min_contig_len_;
}

void BWAISCounter::ProcessAlignments(const MappedPositionT& l, const MappedPositionT& r) {
    ++mapped_count_;

    int is = r.pos - l.pos;
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

// ==== pair info index filler overloads ====
EdgePair BWAIndexFiller::ConjugatePair(EdgePair ep) const {
    return make_pair(g_.conjugate(ep.second), g_.conjugate(ep.first));
}

void BWAIndexFiller::ProcessAlignments(const MappedPositionT& l, const MappedPositionT& r) {
    EdgePair ep{l.e, r.e};
    TRACE("Lpos " << l.pos << ", Rpos " << r.pos);
    int edge_distance = (int) lib_.data().mean_insert_size  - r.pos + l.pos;
    TRACE("Distance " << edge_distance);

    paired_index_.Add(ep.first, ep.second, omnigraph::de::RawPoint(edge_distance, 1.0));
}

bool BWAIndexFiller::CheckAlignments(const MappedPositionT& l, const MappedPositionT& r) {
    return g_.length(l.e) >= min_contig_len_ && g_.length(r.e) >= min_contig_len_;
}


//Main class realization
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
    index_line = path::screen_whitespaces(index_line);
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
    run_command = path::screen_whitespaces(run_command);
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

    //Left and right reads are stored in maps until pair is detected
    unordered_map<string, MapperReadT> left_reads;
    unordered_map<string, MapperReadT> right_reads;
    size_t counter = 0;
    //Check for duplicating read IDs
    bool left_duplicated = false;
    bool right_duplicated = false;

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
            l_name = left_read.name();
            if (left_read.is_properly_aligned()) {
                TRACE("Left read " << l_name);
                left_data = MapperReadT(string(lf.get_contig_name(left_read.contig_id())),
                                        left_read.pos(),
                                        left_read.data_len(),
                                        left_read.strand(),
                                        left_read.cigar());
            }
            else if (!left_read.is_main_alignment()) {
                //If not primary alignment ignore mapping
                TRACE("Ignoring left read");
                l_name = "";
            }
        }
        if (!rf.eof()) {
            rf >> right_read;
            r_name = right_read.name();
            if (right_read.is_properly_aligned()) {
                TRACE("Right read " << r_name);
                right_data = MapperReadT(string(rf.get_contig_name(right_read.contig_id())),
                                         right_read.pos(),
                                         right_read.data_len(),
                                         right_read.strand(),
                                         right_read.cigar());
            }
            else if (!right_read.is_main_alignment()) {
                //If not primary alignment ignore mapping
                TRACE("Ignoring right read");
                r_name = "";
            }
        }

        //Think about custom read names
        if (l_name == r_name) {
            TRACE("Equal processing");
            //Process immideately if ids are equal in both SAM entries
            processor.ProcessPairedRead(left_data, right_data);
            VERBOSE_POWER2(++counter, "Processed " << counter << " paired reads");
            continue;
        }

        if (r_name != "") {
            auto it = left_reads.find(r_name);
            if (it != left_reads.end())  {
                //Right read's mate found in map
                TRACE("Right read's mate found, processing");
                processor.ProcessPairedRead(it->second, right_data);
                VERBOSE_POWER2(++counter, "Processed " << counter << " paired reads");
                //Remove mate as used
                left_reads.erase(it);
            }
            else {
                TRACE("Right read's mate not found, adding to map");
                if (right_reads.count(r_name) == 0) {
                    //Insert read without mate for further analysis
                    //TODO inspect map size and performance
                    right_reads.emplace(r_name, right_data);
                } else {
                    DEBUG("Right read " << r_name << " is duplicated!");
                    //Report duplication
                    right_duplicated = true;
                }
            }
        }

        if (l_name != "") {
            auto it = right_reads.find(l_name);
            if (it != right_reads.end()) {
                //Left read's mate found in map
                TRACE("Left read's mate found, processing");
                processor.ProcessPairedRead(left_data, it->second);
                VERBOSE_POWER2(++counter, "Processed " << counter << " paired reads");
                //Remove mate as used
                right_reads.erase(it);
            }
            else {
                TRACE("Left read's mate not found, adding to map");
                if (left_reads.count(l_name) == 0) {
                    //Insert read without mate for further analysis
                    //TODO inspect map size and performance
                    left_reads.emplace(l_name, left_data);
                } else {
                    DEBUG("Left read " << r_name << " is duplicated!");
                    //Report duplication
                    left_duplicated = true;
                }

            }
        }
    }

    if (left_duplicated)
        WARN("SAM file " << left_sam << " contains duplicated read ids");
    if (right_duplicated)
        WARN("SAM file " << right_sam << " contains duplicated read ids");
}

bool BWAPairInfoFiller::Init() {
    if (!index_constructed_) {
        INFO("Initializing bwa pair info counter, working dir " << work_dir_);
        path::make_dir(base_dir_);
        work_dir_ = path::make_temp_dir(base_dir_, "");
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
    //Initialize if needed
    Init();
    string lib_dir =  path::append_path(work_dir_, ToString(lib_index));
    path::make_dir(lib_dir);
    vector<pair<string, string>> sam_files;
    bool result = false;

    INFO("Mapping lib #" << lib_index << " using BWA");
    if (!AlignLib(lib, path::append_path(lib_dir, "single"), sam_files)) {
        WARN("Failed to align lib #" << lib_index);
        return false;
    }

    INFO("Estimating insert size for library #" << lib_index);
    BWAISCounter counter(lib, edge_id_map_, g_, counter_edge_len);
    for (const auto& sam_pair : sam_files) {
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
        paired_index.Init();

        BWAIndexFiller filler(lib, edge_id_map_, g_, paired_index, index_filler_edge_len);
        for (const auto& sam_pair : sam_files) {
            ProcessSAMFiles(sam_pair.first, sam_pair.second, filler);
        }
        result = true;
    }
    if (remove_tmp_files_)
        path::remove_dir(lib_dir);
    return result;
}


}
