#include "dataset_processor.hpp"
#include <iostream>
#include "io/ireader.hpp"
#include "io/converting_reader_wrapper.hpp"
#include "io/file_reader.hpp"
#include "path_helper.hpp"
#include <omp.h>
using namespace std;
namespace corrector{


void DatasetProcessor::SplitGenome(const string &genome, const string &genome_splitted_dir){
    auto cont_map = GetContigs(genome);
    for (auto &p :cont_map) {
    	string contig_name = p.first;
    	string contig_seq = p.second;
    	string full_path = genome_splitted_dir + "/" + contig_name+ ".fasta";
    	string out_full_path = 	 full_path.substr(0, full_path.length() - 5) + "ref.fasta";
    	string sam_filename =  full_path.substr(0, full_path.length() - 5) + "pair.sam";
    	all_contigs[contig_name].input_contig_filename = full_path;
    	all_contigs[contig_name].output_contig_filename = out_full_path;
    	all_contigs[contig_name].sam_filename = sam_filename;
    	all_contigs[contig_name].contig_length = contig_seq.length();
    	all_contigs[contig_name].buffered_reads.clear();
		PutContig(full_path, contig_name, contig_seq);
    	DEBUG("full_path "+full_path)
    }
}

//contigs - set of aligned contig names
void DatasetProcessor::GetAlignedContigs(const string &read, set<string> &contigs) const {
	vector<string> arr = split(read, '\t');
	if (arr.size() > 5) {
		if (arr[2] != "*" && stoi(arr[4]) > 0) {
// here can be multuple aligned parsing if neeeded;
			contigs.insert(arr[2]);
		}
	}

}


void DatasetProcessor::SplitSingleLibrary(const string &all_reads_filename, const size_t lib_count){
	ifstream fs(all_reads_filename);
	while (! fs.eof()){
		set<string> contigs;
		string r1;
		getline(fs, r1);
		if (r1[0] =='@') continue;
		GetAlignedContigs(r1, contigs);
		for (string contig: contigs) {
			VERIFY_MSG(all_contigs.find(contig) != all_contigs.end(), "wrong contig name in SAM file header: " + contig);
			BufferedOutputRead(r1, contig, lib_count);

		}
	}
	FlushAll(lib_count);
}

void DatasetProcessor::FlushAll(const size_t lib_count)  {
	for (auto &ac :all_contigs) {
		auto stream = new ofstream(ac.second.sam_filenames[lib_count].first.c_str(), std::ios_base::app | std::ios_base::out);
		for (string read: ac.second.buffered_reads){
			*stream << read;
			*stream << '\n';
		}
		stream->close();
		ac.second.buffered_reads.clear();
	}
}

void DatasetProcessor::BufferedOutputRead(const string &read, const string &contig_name, const size_t lib_count) {
	all_contigs[contig_name].buffered_reads.push_back(read);
	buffered_count ++;
	if (buffered_count % buff_size == 0) {
		if (buffered_count % (10*buff_size) == 0)
			INFO("processed " << buffered_count << "reads, flushing");
		FlushAll(lib_count);
	}
}

void DatasetProcessor::SplitPairedLibrary(const string &all_reads_filename, const size_t lib_count){
	ifstream fs(all_reads_filename);
	while (! fs.eof()){
		set<string> contigs;
		string r1;
		string r2;
		getline(fs, r1);
		if (r1[0] =='@') continue;
		getline(fs, r2);
		GetAlignedContigs(r1, contigs);
		GetAlignedContigs(r2, contigs);
		for (string contig: contigs) {
		//	VERIFY_MSG(all_contigs.find(contig) != all_contigs.end(), "wrong contig name in SAM file header: " + contig);
			if (all_contigs.find(contig)!= all_contigs.end()) {
				BufferedOutputRead(r1, contig, lib_count);
				BufferedOutputRead(r2, contig, lib_count);
			}
		}
	}
	FlushAll(lib_count);
}

string DatasetProcessor::GetLibDir(const size_t lib_count) const{
	return corr_cfg::get().work_dir + "/lib" + to_string(lib_count);
}

string DatasetProcessor::RunPairedBwa(const string &left, const string &right, const size_t lib) const{
	string slib = to_string(lib);
	string cur_dir = corr_cfg::get().work_dir + "/lib" + slib;
	string cur_line = "mkdir " + cur_dir;
	int run_res = system(cur_line.c_str());
	string tmp1_sai_filename = cur_dir +"/tmp1.sai";
	string tmp2_sai_filename = cur_dir +"/tmp2.sai";
	string tmp_sam_filename = cur_dir + "/tmp.sam";
	string isize_txt_filename = cur_dir +"/isize.txt";
	string tmp_file = cur_dir + "/bwa.flood";
	string index_line = corr_cfg::get().bwa + " index " + "-a " + "is "+ genome_file + " " + " 2>" + tmp_file;;
	INFO("Running bwa index ...: " << index_line);
	run_res = system(index_line.c_str());
	if (run_res != 0) {
		INFO("bwa failed, skipping sublib");
		return "";
	}
	string nthreads = to_string(corr_cfg::get().max_nthreads);
	string left_line =  corr_cfg::get().bwa + " aln " + genome_file +" "+ left + " -t " + nthreads +
	 				   " -O 7 -E 2 -k 3 -n 0.08 -q 15 > " +tmp1_sai_filename + " 2>" + tmp_file;

	string right_line =  corr_cfg::get().bwa + " aln " + genome_file +" "+ right + " -t " + nthreads +
	 				   " -O 7 -E 2 -k 3 -n 0.08 -q 15 > " +tmp2_sai_filename + " 2>" + tmp_file;
	INFO("Running bwa aln ...: " << left_line);
	run_res += system(left_line.c_str());
	INFO("Running bwa aln ...: " <<right_line);
	run_res += system(right_line.c_str());
	if (run_res != 0) {
		INFO("bwa failed, skipping sublib");
		return "";
	}

	string last_line =  corr_cfg::get().bwa + " sampe "+ genome_file + " " + tmp1_sai_filename + " " +
	 				   tmp2_sai_filename + " " + left + " " + right + "  > " +  tmp_sam_filename + " 2>" + isize_txt_filename;
	INFO("Running bwa sampe ...:" << last_line);
	run_res = system(last_line.c_str());
	if (run_res != 0) {
		INFO("bwa failed, skipping sublib");
		return "";
	}
	return tmp_sam_filename;
}

string DatasetProcessor::RunSingleBwa(const string &single,const size_t lib) const{
	int run_res = 0;
	string slib = to_string(lib);
	string cur_dir = corr_cfg::get().work_dir + "/lib" + slib;
	string cur_line = "mkdir " + cur_dir;
	run_res = system(cur_line.c_str());
	string tmp_sai_filename = cur_dir +"/tmp.sai";

	string tmp_sam_filename = cur_dir + "/tmp.sam";
	string isize_txt_filename = cur_dir +"/isize.txt";
	string tmp_file = cur_dir + "/bwa.flood";
	string index_line = corr_cfg::get().bwa + " index " + "-a " + "is "+ genome_file + " 2 "+ "2>" + tmp_file;
	INFO("Running bwa index ...: " << index_line);
	run_res= system(index_line.c_str());
	if (run_res != 0) {
		INFO("bwa failed, skipping sublib");
		return "";
	}
	string nthreads = to_string(corr_cfg::get().max_nthreads);
	string single_sai_line =  corr_cfg::get().bwa + " aln " + genome_file +" "+ single+ " -t " + nthreads +
		 				   " -O 7 -E 2 -k 3 -n 0.08 -q 15 > " +tmp_sai_filename + " 2>" + tmp_file;

	INFO("Running bwa aln ...:" + single_sai_line);

	run_res = system(single_sai_line.c_str());
	if (run_res != 0) {
		INFO("bwa failed, skipping sublib");
		return "";
	}
	string last_line =  corr_cfg::get().bwa + " samse "+ genome_file + " " + tmp_sai_filename + " " +
	 				    single +"  > " +  tmp_sam_filename + " 2>" + isize_txt_filename;
	INFO("Running bwa samse ...:" << last_line);
	run_res = system(last_line.c_str());
	if (run_res != 0) {
		INFO("bwa failed, skipping sublib");
		return "";
	}
	return tmp_sam_filename;

}


void DatasetProcessor::PrepareContigDirs(const size_t lib_count) {

	string out_dir = GetLibDir(lib_count);
	for(auto &ac: all_contigs ){
		auto contig_name = ac.first;
		string header = "@SQ\tSN:"+ contig_name+"\tLN:"+ to_string(all_contigs[contig_name].contig_length);
		BufferedOutputRead(header, contig_name, lib_count);
		string out_name = out_dir + "/" + ac.first + ".sam";
		ac.second.sam_filenames.push_back(make_pair(out_name, unsplitted_sam_files[lib_count].second));
	}
	FlushAll(lib_count);
}
void DatasetProcessor::ProcessDataset() {
	size_t lib_num = 0;
	INFO("Splitting genome...");
	INFO("Genome_file: " + genome_file);
	SplitGenome(genome_file, work_dir);
    for (size_t i = 0 ; i < corr_cfg::get().dataset.lib_count(); ++i) {
    	if (corr_cfg::get().dataset[i].type() == io::LibraryType::PairedEnd || corr_cfg::get().dataset[i].type() == io::LibraryType::HQMatePairs || corr_cfg::get().dataset[i].type() == io::LibraryType::SingleReads) {
    		for (auto iter = corr_cfg::get().dataset[i].paired_begin(); iter !=  corr_cfg::get().dataset[i].paired_end(); iter++){
    			INFO("Processing paired sublib of number " << lib_num);
    			string left = iter->first;
    			string right = iter->second;
    			INFO(left + " " + right);
    			string samf = RunPairedBwa(left, right, lib_num);
    			//INFO(RunPairedBwa(left, right, lib_num));
    			if (samf !="") {
        			INFO("Adding samfile " << samf);

    				unsplitted_sam_files.push_back(make_pair(samf, "paired"));
    				PrepareContigDirs(lib_num);
    				SplitPairedLibrary(samf, lib_num);
            		lib_num++;
    			} else {
    				WARN("Failed to align paired reads " << left << " and " << right);
    			}
    		}
    		for (auto iter = corr_cfg::get().dataset[i].single_begin(); iter !=  corr_cfg::get().dataset[i].single_end(); iter++){
    			INFO("Processing single sublib of number " << lib_num);
    			string left = *iter;
				INFO(left);
				string samf = RunSingleBwa(left, lib_num);
				if (samf != "") {
					INFO("Adding samfile "<< samf);
					unsplitted_sam_files.push_back(make_pair(samf, "single"));
					PrepareContigDirs(lib_num);
					SplitSingleLibrary(samf, lib_num);
					lib_num++;
				} else {
					WARN("Failed to align single reads " << left );
				}
			}
    	}
    }
	INFO("Processing contigs");
	vector<pair<size_t, string> > ordered_contigs;
	ContigInfoMap all_contigs_copy;
	for (auto ac: all_contigs) {
		ordered_contigs.push_back(make_pair(ac.second.contig_length, ac.first));
		all_contigs_copy.insert(ac);
	}
	size_t cont_num = ordered_contigs.size();
	sort(ordered_contigs.begin(), ordered_contigs.end());
	reverse(ordered_contigs.begin(), ordered_contigs.end());
# pragma omp parallel for shared( all_contigs_copy, ordered_contigs) num_threads(nthreads) schedule(dynamic,1)
	for (size_t i = 0; i < cont_num; i++ ) {

		if ( ordered_contigs[i].first > 20000) {
			INFO("Processing contig" << ordered_contigs[i].second << " of length " << ordered_contigs[i].first <<" in thread number " << omp_get_thread_num());
		}
		ContigProcessor pc(all_contigs_copy[ordered_contigs[i].second].sam_filenames, all_contigs_copy[ordered_contigs[i].second].input_contig_filename);

		pc.ProcessMultipleSamFiles();
	}
	INFO("Gluing processed contigs");
	GlueSplittedContigs(output_contig_file);
}

void DatasetProcessor::GlueSplittedContigs(string &out_contigs_filename){
	ofstream of_c(out_contigs_filename, std::ios_base::binary);
	for (auto ac : all_contigs) {
		string a = ac.second.output_contig_filename;
		ifstream a_f(a, std::ios_base::binary);
		of_c << a_f.rdbuf();
	}
}



};
