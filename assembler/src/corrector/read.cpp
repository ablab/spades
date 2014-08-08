#include "read.hpp"
#include "variants_table.hpp"

#include "logger/log_writers.hpp"

using namespace std;


namespace corrector {

int SingleSamRead::CountPositions(unordered_map<size_t, position_description> &ps, const size_t &contig_length) const {

    if (get_contig_id() < 0) {
        DEBUG("not this contig");
        return -1;
    }
    if (data_->core.qual == 0) {
        DEBUG("zero qual");
        return -1;
    }
    int pos = data_->core.pos;
    if (pos < 0) {
        WARN("Negative position " << pos << " found on read " << get_name() << ", skipping");
        return -1;
    }
    size_t position = size_t(pos);
    int mate = 1;  // bonus for mate mapped can be here;
    size_t l_read = get_data_len();
    size_t l_cigar = get_cigar_len();

    int aligned_length = 0;
    uint32_t *cigar = bam1_cigar(data_);
    //* in cigar;
    if (l_cigar == 0)
        return -1;
    if (bam_cigar_opchr(cigar[0]) == '*')
        return -1;
    for (size_t i = 0; i < l_cigar; i++)
        if (bam_cigar_opchr(cigar[i]) == 'M')
            aligned_length += bam_cigar_oplen(cigar[i]);
//It's about bad aligned reads, but whether it is necessary?
    double read_len_double = (double) l_read;
    if ((aligned_length < min(read_len_double * 0.4, 40.0)) && (position > read_len_double / 2) && (contig_length > read_len_double / 2 + (double) position)) {
        return -1;
    }
    int state_pos = 0;
    int shift = 0;
    size_t skipped = 0;
    size_t deleted = 0;
    string insertion_string = "";
    auto seq = bam1_seq(data_);
    for (size_t i = 0; i < l_read; i++) {
        DEBUG(i << " " << position << " " << skipped);
        if (shift + bam_cigar_oplen(cigar[state_pos]) <= i) {
            shift += bam_cigar_oplen(cigar[state_pos]);
            state_pos += 1;
        }
        if (insertion_string != "" and bam_cigar_opchr(cigar[state_pos]) != 'I') {
            VERIFY(i + position >= skipped + 1);
            size_t ind = i + position - skipped - 1;
            if (ind >= contig_length)
                break;
            ps[ind].insertions[insertion_string] += 1;
            insertion_string = "";
        }
        char cur_state = bam_cigar_opchr(cigar[state_pos]);
        if (cur_state == 'M') {
            VERIFY(i >= deleted);
            if (i + position < skipped) {
                WARN(i << " " << position << " " << skipped);
                INFO(get_name());
            }
            VERIFY(i + position >= skipped);

            size_t ind = i + position - skipped;
            size_t cur = var_to_pos[(int) bam_nt16_rev_table[bam1_seqi(seq, i - deleted)]];
            if (ind >= contig_length)
                continue;
            ps[ind].votes[cur] = ps[ind].votes[cur] + mate;

        } else {
            if (cur_state == 'I' || cur_state == 'H' || cur_state == 'S' ) {
                if (cur_state == 'I') {
                    if (insertion_string == "") {
                        size_t ind = i + position - skipped - 1;
                        if (ind >= contig_length)
                            break;
                        ps[ind].votes[Variants::Insertion] += mate;
                    }
                    insertion_string += bam_nt16_rev_table[bam1_seqi(seq, i - deleted)];
                }
                skipped += 1;
            } else if (bam_cigar_opchr(cigar[state_pos]) == 'D') {
                if (i + position - skipped >= contig_length)
                    break;
                ps[i + position - skipped].votes[Variants::Deletion] += mate;
                deleted += 1;
            }
        }
    }
    if (insertion_string != "" and bam_cigar_opchr(cigar[state_pos]) != 'I') {
        VERIFY(l_read + position >= skipped + 1);
        size_t ind = l_read + position - skipped - 1;
        if (ind < contig_length) {
            ps[ind].insertions[insertion_string] += 1;
        }
        insertion_string = "";
    }
    return 0;
}

string SingleSamRead::get_cigar() const {
    uint32_t *cigar = bam1_cigar(data_);
    string res;
    res.reserve(data_->core.n_cigar);
    for (size_t k = 0; k < data_->core.n_cigar; ++k) {
        res += std::to_string(bam_cigar_oplen(cigar[k]));
        res += bam_cigar_opchr(cigar[k]);

    }
    return res;
}

string SingleSamRead::get_name() const {
    string res(bam1_qname(data_));
    return res;
}

string SingleSamRead::get_seq() const {
    string res = "";
    auto b = bam1_seq(data_);
    for (int k = 0; k < data_->core.l_qseq; ++k) {
        res += bam_nt16_rev_table[bam1_seqi(b, k)];
    }
    return res;
}

int PairedSamRead::CountPositions(unordered_map<size_t, position_description> &ps, const size_t &contig_length) const {

    TRACE("starting pairing");
    int t1 = r1.CountPositions(ps, contig_length);
    unordered_map<size_t, position_description> tmp;
    int t2 = r2.CountPositions(tmp, contig_length);
    //overlaps.. multimap? Look on qual?
    if (ps.size() == 0 || tmp.size() == 0) {
        //We do not need paired reads which are not really paired
        ps.clear();
        return -1;
    }
    TRACE("counted, uniting maps of " << tmp.size() << " and " << ps.size());
    ps.insert(tmp.begin(), tmp.end());
    TRACE("united");
    return t1 + t2;
}

}
;
