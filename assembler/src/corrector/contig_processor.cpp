#include "contig_processor.hpp"
#include "include.hpp"
#include "io/ireader.hpp"
#include "io/converting_reader_wrapper.hpp"
#include "io/file_reader.hpp"

namespace corrector{

void ContigProcessor::read_contig() {
	auto res= GetContigs(contig_file);
	VERIFY_MSG(res.size()== 1, "Not one contig in file with unique contigs");

	contig = res.begin()->second;
	contig_name = res.begin()->first;
	if (debug_info) {
		INFO("Processing contig of length " << contig.length());
	}
	DEBUG("Processing contig of length " << contig.length());
//extention is always "fasta"
	output_contig_file = contig_file.substr(0, contig_file.length() - 5) + "ref.fasta";
	charts.resize(contig.length());
}

void ContigProcessor::UpdateOneRead(const SingleSamRead &tmp, MappedSamStream &sm){
	unordered_map<size_t, position_description> all_positions;
	if (tmp.get_contig_id() < 0) {
		return;
	}
	string cur_s = sm.get_contig_name(tmp.get_contig_id());
	if (cur_s != contig_name) {
		return;
	}
	int error_num = tmp.CountPositions(all_positions, contig);
	if (error_num > 19) error_counts[20] ++;
	else if (error_num >=0) error_counts[error_num] ++;
	for (auto &pos :all_positions) {
		if ((int)pos.first >=0 && pos.first < contig.length()) {
			charts[pos.first].update(pos.second);
		}
	}
}

//returns: number of changed nucleotides;
int ContigProcessor::UpdateOneBase(size_t i, stringstream &ss, const unordered_map<size_t, position_description> &interesting_positions){
	char old = (char) toupper(contig[i]);
	size_t maxi = charts[i].FoundOptimal(contig[i]);
	if (interesting_positions.find(i) != interesting_positions.end()) {
		size_t maxj = interesting_positions.find(i)->second.FoundOptimal(contig[i]);
		if (maxj != maxi) {
			DEBUG("Interesting positions differ with majority!");
			DEBUG("On position " << i << "  old: " << old <<" majority: "<<pos_to_var[maxi] << "interesting: " << pos_to_var[maxj]);
			if (corr_cfg::get().strategy != "majority_only")
				maxi = maxj;
		}
	}
	if (old != pos_to_var[maxi]) {
		DEBUG("On position " << i << " changing " << old <<" to "<<pos_to_var[maxi]);
		DEBUG(charts[i].str());
		if (maxi < DELETION) {
			ss <<pos_to_var[maxi];
			return 1;
		} else if (maxi == DELETION) {
			return 1;
		} else if (maxi == INSERTION) {
			string maxj = "";
			//first base before insertion;
			size_t new_maxi = var_to_pos[(int)contig[i]];
			int new_maxx = charts[i].votes[new_maxi];
			for (size_t k = 0; k < MAX_VARIANTS; k++) {
				if (new_maxx < charts[i].votes[k] && (k != INSERTION) && (k != DELETION)) {
					new_maxx = charts[i].votes[k];
					new_maxi = k;
				}
			}
			ss <<pos_to_var[new_maxi];
			int max_ins = 0;
			for (auto &ic :charts[i].insertions) {
				if (ic.second > max_ins){
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
			//something starnge happened
			WARN("While processing base " << i << " unknown decision was made");
			return 0;
		}
	} else {
		ss << old;
		return 0;
	}
}

void ContigProcessor::process_multiple_sam_files() {
	error_counts.resize(20);
	DEBUG("working with " << sam_files.size() << " sublibs");
	for (auto &sf : sam_files){
		MappedSamStream sm (sf.first);
		while (!sm.eof()) {
			SingleSamRead tmp;
			sm >> tmp;
			UpdateOneRead(tmp, sm);
		}
	}
	stringstream err_str;
	for(int i = 0; i < 20; i ++)
		err_str << error_counts[i] << " ";
	size_t interesting = ipp.FillInterestingPositions(charts);
	if (debug_info) {
		INFO("Error counts for contig "<<contig_name<< " in thread " << omp_get_thread_num() <<" : " << err_str.str() );
	} else {
		DEBUG("Error counts for contig "<<contig_name<<" : " << err_str.str() );
	}
	DEBUG("Interesting size: " << interesting);
	for (auto &sf : sam_files){
		MappedSamStream sm (sf.first);
		while (!sm.eof()) {
			unordered_map<size_t, position_description> ps;
			if (sf.second == "paired") {
				PairedSamRead tmp;
				sm >>tmp;
				tmp.CountPositions(ps, contig);
			} else {
				SingleSamRead tmp;
				sm >>tmp;
				tmp.CountPositions(ps, contig);
			}
			TRACE("updating interesting read..");
			ipp.UpdateInterestingRead(ps);
		}
	}
	ipp.UpdateInterestingPositions();
	unordered_map<size_t, position_description> interesting_positions = ipp.get_weights();
	stringstream s_new_contig;
	int total_changes = 0;
	for (size_t i = 0; i < contig.length(); i ++) {
		DEBUG(charts[i].str());
		total_changes += UpdateOneBase(i, s_new_contig, interesting_positions);
	}
	if (debug_info || total_changes * 5 > (int) contig.length()) {
		INFO("Contig " << contig_name <<" processed with " << total_changes << " changes");
	} else {
		DEBUG("Contig " << contig_name <<" processed with " << total_changes << " changes");
	}
	contig_name = ContigRenameWithLength(contig_name, s_new_contig.str().length());
	PutContig(output_contig_file, contig_name, s_new_contig.str());
}

};
