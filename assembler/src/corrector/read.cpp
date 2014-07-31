#include "read.hpp"
#include "include.hpp"
using namespace std;

namespace corrector {

void position_description::update(const position_description &another){
	for (size_t i = 0; i < MAX_VARIANTS; i++)
		votes[i] += another.votes[i];
	for (auto &ins : another.insertions)
		insertions[ins.first] += ins.second;
}
string position_description::str() const{
	stringstream ss;
	for (int i = 0; i < MAX_VARIANTS; i++ ){
		ss << pos_to_var[i];
		ss <<  ": " << votes[i]<<"; ";

	}
	return ss.str();
}
size_t position_description::FoundOptimal(char current) const{
	size_t maxi = var_to_pos[(size_t) current];
	int maxx = votes[maxi];
	for (size_t j = 0; j < MAX_VARIANTS; j++) {
		//1.5 because insertion goes _after_ match
		//std::min_element
		if (maxx < votes[j] || (j == INSERTION && maxx * 2 <votes[j] * 3)) {
			maxx = votes[j];
			maxi = j;
		}
	}
	return maxi;
}
void position_description::clear() {
	for (size_t i = 0; i < MAX_VARIANTS; i++) {
		votes[i] = 0;
		insertions.clear();
	}
}

size_t SingleSamRead::DataLen() const{
	return data_.core.l_qseq;
}
size_t SingleSamRead::CigarLen() const{
	return data_.core.n_cigar;
}
int SingleSamRead::get_contig_id() const{
	return data_.core.tid;
}
void SingleSamRead::set_data(bam1_t *seq_) {
	//TODO: delete
	bam1_t *new_seq = bam_dup1(seq_);
	data_ = *new_seq;
}
int SingleSamRead::CountPositions(unordered_map <size_t, position_description> &ps, const string &contig) const{
	size_t contig_length = contig.length();
	int error_num = 0;
	if (get_contig_id() < 0) {
		DEBUG("not this contig");
		return -1;
	}
	if (data_.core.qual == 0) {
		DEBUG("zero qual");
		return -1;
	}
	int pos = data_.core.pos;
	if (pos < 0) {
		WARN("Negative position " << pos << " found on read " << GetName() <<", skipping");

		return -1;
	}
	size_t position = size_t(pos);
	int mate = 1; // bonus for mate mapped can be here;
	size_t l_read = DataLen();
	size_t l_cigar = CigarLen();

	set<char> to_skip = {'S', 'I', 'H'};
	int  aligned_length = 0;
	uint32_t *cigar = bam1_cigar(&data_);
   //* in cigar;
	if (l_cigar == 0)
		return -1;
	if (bam_cigar_opchr(cigar[0]) =='*')
		return -1 ;
	for (size_t i = 0; i <l_cigar; i++)
		if (bam_cigar_opchr(cigar[i]) =='M')
			aligned_length +=  bam_cigar_oplen(cigar[i]);
//It's about bad aligned reads, but whether it is necessary?
	double read_len_double = (double) l_read;
	if ((aligned_length < min(read_len_double* 0.4, 40.0)) && (position > read_len_double/ 2) && (contig_length  > read_len_double/ 2 + (double) position) ) {
		return -1;
	}
	int state_pos = 0;
	int shift = 0;
	size_t skipped = 0;
	size_t deleted = 0;
	string insertion_string = "";

	auto seq = bam1_seq(&data_);
	for (size_t i = 0; i < l_read; i++) {
		DEBUG(i << " " << position << " " << skipped);
		if (shift +  bam_cigar_oplen(cigar[state_pos]) <= i){
			shift +=  bam_cigar_oplen(cigar[state_pos]);
			state_pos += 1;
		}
		if (insertion_string != "" and bam_cigar_opchr(cigar[state_pos]) != 'I'){
			VERIFY(i + position >= skipped + 1);
			size_t ind = i + position - skipped - 1;
			if (ind >= contig_length)
				break;
			ps[ind].insertions[insertion_string] += 1;
			insertion_string = "";
		}
		char cur_state = bam_cigar_opchr(cigar[state_pos]);
		if (cur_state == 'M'){

			VERIFY(i >= deleted);
			if (i + position < skipped) {
				WARN(i << " " << position <<" "<< skipped);
				INFO(GetName());
			}
			VERIFY (i + position >= skipped);

			size_t ind = i + position - skipped;
			int cur = var_to_pos[(int)bam_nt16_rev_table[bam1_seqi(seq, i - deleted)]];
			if (ind >= contig_length)
				continue;
			ps[ind].votes[cur] = ps[ind].votes[cur] +  mate;//t_mate
			if (contig[ind] != pos_to_var[cur]) error_num ++;
		} else {
			if (to_skip.find(cur_state) != to_skip.end()){
				if (cur_state == 'I'){
					if (insertion_string == "") {
						size_t ind = i + position - skipped - 1;
						if (ind >= contig_length)
							break;
						ps[i + position - skipped - 1].votes[INSERTION] += mate;
					}
					insertion_string += bam_nt16_rev_table[bam1_seqi(seq, i - deleted)];
					error_num ++;
				}
				skipped += 1;
		   } else if (bam_cigar_opchr(cigar[state_pos]) == 'D') {
			   if (i + position - skipped >= contig_length)
				   break;
				ps[i + position - skipped].votes[DELETION] += mate;
				deleted += 1;
				error_num ++;
		   }
		}
	}
	if (insertion_string != "" and bam_cigar_opchr(cigar[state_pos]) != 'I'){
		VERIFY(l_read + position >= skipped + 1);
		size_t ind = l_read + position - skipped - 1;
		if (ind < contig_length) {
			ps[ind].insertions[insertion_string] += 1;
		}
		insertion_string = "";
	}
	if (false) {
		INFO("strange read");
		INFO(GetName());
		INFO(GetSeq());
		INFO(GetCigar());

		INFO("position of selected " << position);
		INFO ("first char " << bam_cigar_opchr(cigar[0]));
	}
	return error_num;
}

string SingleSamRead::GetCigar() const {
	uint32_t *cigar = bam1_cigar(&data_);
	string res;
	res.reserve(data_.core.n_cigar);
	for (size_t k = 0; k < data_.core.n_cigar; ++k) {
		res += std::to_string( bam_cigar_oplen(cigar[k]));
		res += bam_cigar_opchr(cigar[k]);

	}
	return res;
}

string SingleSamRead::GetQual() const{
	uint8_t *qual = bam1_qual(&data_);
	for (int i = 0; i < data_.core.l_qseq; ++i) {
		qual[i] = uint8_t (qual[i] + 33);
	}
	string res(reinterpret_cast<const char*>(qual));
	return res;
}

string SingleSamRead::GetName() const{
	string res(bam1_qname(&data_));
	return res;
}

string SingleSamRead::GetSeq() const{
	string res = "";
	auto b = bam1_seq(&data_);
	for (int k = 0; k < data_.core.l_qseq; ++k) {
		res += bam_nt16_rev_table[bam1_seqi(b, k)];
	}
	return res;
}


void PairedSamRead::pair(SingleSamRead &a1, SingleSamRead &a2) {
	r1 = a1; r2 = a2;
}

int PairedSamRead::CountPositions(unordered_map <size_t, position_description> &ps, const string &contig) const{

	TRACE("starting pairing");
	int t1 = r1.CountPositions(ps, contig);
	unordered_map <size_t, position_description> tmp;
	int t2 = r2.CountPositions(tmp, contig);
	//overlaps.. multimap? Look on qual?
	if (ps.size() ==0 || tmp.size() == 0) {
		//We do not need paired reads which are not really paired
		ps.clear();
		return -1 ;
	}
	TRACE("counted, uniting maps of " << tmp.size () << " and " << ps.size());
	ps.insert( tmp.begin(), tmp.end());
	TRACE("united");
	return t1 + t2;
}



};
