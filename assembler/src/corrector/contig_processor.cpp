#include "contig_processor.hpp"
#include "include.hpp"
#include "io/ireader.hpp"
#include "io/converting_reader_wrapper.hpp"
#include "io/file_reader.hpp"

namespace corrector{

void ContigProcessor::read_contig() {
	io::FileReadStream contig_stream(contig_file);
	io::SingleRead ctg;
	contig_stream >> ctg;
	contig = ctg.sequence().str();
	contig_name = ctg.name();
	INFO("Processing contig of length " << contig.length());
//extention is always "fasta"
	output_contig_file = contig_file.substr(0, contig_file.length() - 5) + "ref.fasta";
	charts.resize(contig.length());
}

void ContigProcessor::UpdateOneRead(SingleSamRead &tmp){
	unordered_map<size_t, position_description> all_positions;
	//INFO(tmp.GetName());
	//INFO(tmp.get_contig_id());
	if (tmp.get_contig_id() < 0) {
		return;
	}
	string cur_s = sm.get_contig_name(tmp.get_contig_id());

	if (cur_s != contig_name) {
		WARN("wrong string");
		return;
	}
	int error_num = tmp.CountPositions(all_positions, contig);
	if (error_num > 19) error_counts[20] ++;
	else if (error_num >=0) error_counts[error_num] ++;
	for (auto iter = all_positions.begin(); iter != all_positions.end(); ++iter) {
		if ((int)iter->first >=0 && iter->first < contig.length()) {
			charts[iter->first].update(iter->second);


		}
	}
}
//returns: number of changed nucleotides;
int ContigProcessor::UpdateOneBase(size_t i, stringstream &ss, unordered_map<size_t, position_description> &interesting_positions){
	char old = (char) toupper(contig[i]);
	size_t maxi = charts[i].FoundOptimal(contig[i]);
	if (interesting_positions.find(i) != interesting_positions.end()) {
		size_t maxj = interesting_positions[i].FoundOptimal(contig[i]);
		if (maxj != maxi) {
			INFO("Interesting positions differ with majority!");
			INFO("On position " << i << "  old: " << old <<" majority: "<<pos_to_var[maxi] << "interesting: " << pos_to_var[maxj]);
			maxi = maxj;
		}
	}
	if (old != pos_to_var[maxi]) {
		INFO("On position " << i << " changing " << old <<" to "<<pos_to_var[maxi]);
		INFO(charts[i].str());
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
			for (auto iter = charts[i].insertions.begin(); iter != charts[i].insertions.end(); ++iter) {
				if (iter->second > max_ins){
					max_ins = iter->second;
					maxj = iter->first;
				}
			}
			INFO("most popular insertion: " << maxj);
			ss << maxj;
			if (old == maxj[0]) {
				return (int) maxj.length() - 1;
			} else {
				return (int) maxj.length();
			}
		} else {
			//something starnge happened
			WARN("While processing base " << i << " unknown decision was made");
			return -1;
		}
	} else {
		ss << old;
		return 0;
	}
}

void ContigProcessor::process_sam_file (){
	error_counts.resize(20);
	while (!sm.eof()) {
		SingleSamRead tmp;
		sm >> tmp;
		UpdateOneRead(tmp);
		//	sm >>tmp;
	}
	for(int i = 0; i < 20; i ++)
		cout << error_counts[i] << " ";
	sm.reset();
	size_t interesting = ipp.FillInterestingPositions(charts);
	INFO("interesting size: " << interesting);
	while (!sm.eof()) {
		PairedSamRead tmp;
		unordered_map<size_t, position_description> ps;
		sm >>tmp;
		tmp.CountPositions(ps, contig);
		TRACE("updating interesting read..");
		ipp.UpdateInterestingRead(ps);
	}
	ipp.UpdateInterestingPositions();
	unordered_map<size_t, position_description> interesting_positions = ipp.get_weights();
	stringstream s_new_contig;
	for (size_t i = 0; i < contig.length(); i ++) {
		DEBUG(charts[i].str());
		UpdateOneBase(i, s_new_contig, interesting_positions);
	}
	io::osequencestream oss(output_contig_file);
	contig_name = ContigRenameWithLength(contig_name, s_new_contig.str().length());
	oss << io::SingleRead(contig_name, s_new_contig.str());
}

};
