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
    	//INFO(contig_name);
    	string full_path = genome_splitted_dir + "/" + contig_name+ ".fasta";
    	string out_full_path = 	 full_path.substr(0, full_path.length() - 5) + "ref.fasta";
    	string sam_filename =  full_path.substr(0, full_path.length() - 5) + "pair.sam";
    	//INFO("full_path:" << full_path);
    	all_contigs[contig_name].input_contig_filename = full_path;
    	all_contigs[contig_name].output_contig_filename = out_full_path;
    	all_contigs[contig_name].sam_filename = sam_filename;
    	all_contigs[contig_name].contig_length = contig_seq.length();
    	all_contigs[contig_name].buffered_reads.clear();
		//auto tmp_stream = new ofstream(sam_filename.c_str(), ios::app);
		//tmp_stream->close();
		PutContig(full_path, contig_name, contig_seq);

    	DEBUG("full_path "+full_path)
    	//oss << ctg;

    }
}
//returns names
void DatasetProcessor::GetAlignedContigs(const string &read, set<string> &contigs) const {
	vector<string> arr = split(read, '\t');
	if (arr.size() > 5) {
		if (arr[2] != "*" && stoi(arr[4]) > 0) {
//TODO:: here can be multuple aligned parsing if neeeded;
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

/*
void DatasetProcessor::PrepareWriters(){
//TODO::place for buffered writers;
	for (auto ac :all_contigs){

		all_writers[ac.first] = new ofstream(ac.second.sam_filename.c_str());
	}
}

void DatasetProcessor::CloseWriters(){
//TODO::place for buffered writers;
	for (auto ac :all_contigs){
		all_writers[ac.first]->close();
	}
}


void DatasetProcessor::OutputRead(string &read, string &contig_name) {
	//TODO:: buffered_writer
	*all_writers[contig_name] << read ;
	*all_writers[contig_name] <<  '\n';
}
*/
void DatasetProcessor::FlushAll(const size_t lib_count)  {
	for (auto &ac :all_contigs) {
		auto stream = new ofstream(ac.second.sam_filenames[lib_count].first.c_str(), std::ios_base::app | std::ios_base::out);
		//INFO(ac.second.buffered_reads.size());
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
/*
void DatasetProcessor::SplitHeaders(string &all_reads_filename) {
	ifstream fs(all_reads_filename);
	while (! fs.eof()){
		string r;
		getline(fs, r);
		if (r.length() == 0 || r[0] != '@') {
			break;
		}
		vector<string > arr = split(r, '\t');
		if (arr[0] == "@SQ") {
			VERIFY_MSG(arr.size() > 1, "Invalid .sam header");
			string contig_name = arr[1].substr(3, arr[1].length() - 3) ;
			INFO(contig_name);
			//VERIFY_MSG(all_writers.find(contig_name) != all_writers.end(), "wrong contig name in SAM file header");
			BufferedOutputRead(r, contig_name, lib_count);
		}
	}
}
*/

string DatasetProcessor::RunPairedBwa(const string &left, const string &right, const size_t lib) const{
	string slib = to_string(lib);
	string cur_dir = corr_cfg::get().work_dir + "/lib" + slib;
	string cur_line = "mkdir " + cur_dir;
	system(cur_line.c_str());
	string tmp1_sai_filename = cur_dir +"/tmp1.sai";
	string tmp2_sai_filename = cur_dir +"/tmp2.sai";
	string tmp_sam_filename = cur_dir + "/tmp.sam";
	string isize_txt_filename = cur_dir +"/isize.txt";
	//os.path.join(config["work_dir"], "isize.txt")

	//universal_sys_call([config["bwa"], "index", "-a", "is", config["contigs"], "2"], log)
	string index_line = corr_cfg::get().bwa + " index " + "-a " + "is "+ genome_file + " 2 ";
	INFO("index line: " << index_line);
	system(index_line.c_str());
	string nthreads = to_string(corr_cfg::get().max_nthreads);
	string left_line =  corr_cfg::get().bwa + " aln " + genome_file +" "+ left + " -t " + nthreads +
	 				   " -O 7 -E 2 -k 3 -n 0.08 -q 15 > " +tmp1_sai_filename;

	string right_line =  corr_cfg::get().bwa + " aln " + genome_file +" "+ right + " -t " + nthreads +
	 				   " -O 7 -E 2 -k 3 -n 0.08 -q 15 > " +tmp2_sai_filename;
	INFO(left_line);
	INFO(right_line);
	system(left_line.c_str());
	system(right_line.c_str());
/*	universal_sys_call([config["bwa"], "aln", config["contigs"], config["reads1"], "-t",
				   str(config["t"]), "-O", "7", "-E", "2", "-k", "3", "-n", "0.08", "-q", "15"], log, tmp1_sai_filename)
	universal_sys_call([config["bwa"], "aln", config["contigs"], config["reads2"], "-t",
				   str(config["t"]), "-O", "7", "-E", "2", "-k", "3", "-n", "0.08", "-q", "15"], log, tmp2_sai_filename)

	universal_sys_call([config["bwa"], "sampe", config["contigs"], tmp1_sai_filename,
				   tmp2_sai_filename, config["reads1"], config["reads2"]],
				   None, tmp_sam_filename, isize_txt_filename)*/
	string last_line =  corr_cfg::get().bwa + " sampe "+ genome_file + " " + tmp1_sai_filename + " " +
	 				   tmp2_sai_filename + " " + left + " " + right + "  > " +  tmp_sam_filename + " 2>" + isize_txt_filename;
	INFO(last_line);
	system(last_line.c_str());

	return tmp_sam_filename;
}

string DatasetProcessor::RunSingleBwa(const string &single,const size_t lib) const{
	string slib = to_string(lib);
	string cur_dir = corr_cfg::get().work_dir + "/lib" + slib;
	string cur_line = "mkdir " + cur_dir;
	system(cur_line.c_str());
	string tmp_sai_filename = cur_dir +"/tmp.sai";

	string tmp_sam_filename = cur_dir + "/tmp.sam";
	string isize_txt_filename = cur_dir +"/isize.txt";
		//os.path.join(config["work_dir"], "isize.txt")

		//universal_sys_call([config["bwa"], "index", "-a", "is", config["contigs"], "2"], log)
	string index_line = corr_cfg::get().bwa + " index " + "-a " + "is "+ genome_file + " 2 ";
	INFO("index line: " << index_line);
	system(index_line.c_str());
	string nthreads = to_string(corr_cfg::get().max_nthreads);
	string single_sai_line =  corr_cfg::get().bwa + " aln " + genome_file +" "+ single+ " -t " + nthreads +
		 				   " -O 7 -E 2 -k 3 -n 0.08 -q 15 > " +tmp_sai_filename;

	INFO(single_sai_line);

	system(single_sai_line.c_str());
	string last_line =  corr_cfg::get().bwa + " samse "+ genome_file + " " + tmp_sai_filename + " " +
	 				    single +"  > " +  tmp_sam_filename + " 2>" + isize_txt_filename;
	INFO(last_line);
	system(last_line.c_str());
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
}
void DatasetProcessor::ProcessDataset() {
	size_t lib_num = 0;
	INFO("Splitting genome");
	INFO("genome_file: " + genome_file);
	SplitGenome(genome_file, work_dir);
	//PrepareWriters();
	//SplitHeaders(sam_file);

    for (size_t i = 0 ; i < corr_cfg::get().dataset.lib_count(); ++i) {
    	//INFO(corr_cfg::get().dataset[i].type());
    	if (corr_cfg::get().dataset[i].type() == io::LibraryType::PairedEnd || corr_cfg::get().dataset[i].type() == io::LibraryType::HQMatePairs) {
    		INFO("paired");
    		for (auto iter = corr_cfg::get().dataset[i].paired_begin(); iter !=  corr_cfg::get().dataset[i].paired_end(); iter++){
    			string left = iter->first;
    			string right = iter->second;
    			INFO(left + " " + right);
    			string samf = RunPairedBwa(left, right, lib_num);
    			//INFO(RunPairedBwa(left, right, lib_num));
    			INFO("adding samfile " << samf);
    			unsplitted_sam_files.push_back(make_pair(samf, "paired"));
    			PrepareContigDirs(lib_num);
    			SplitPairedLibrary(samf, lib_num);

        		lib_num++;
    		}
    		for (auto iter = corr_cfg::get().dataset[i].single_begin(); iter !=  corr_cfg::get().dataset[i].single_end(); iter++){
				string left = *iter;
				//string right = iter->second;
				INFO(left);
				string samf = RunSingleBwa(left, lib_num);
				INFO("adding samfile "<< samf);
				unsplitted_sam_files.push_back(make_pair(samf, "single"));
    			PrepareContigDirs(lib_num);
    			SplitSingleLibrary(samf, lib_num);

				lib_num++;

			}
    	} else if (corr_cfg::get().dataset[i].type() == io::LibraryType::SingleReads) {
    		INFO("single");

    		//for (auto iter = corr_cfg::get().dataset[i].paired_begin(); iter !=  corr_cfg::get().dataset[i].single_end(); iter++){
    	}
//
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

# pragma omp parallel for shared(all_contigs_copy, ordered_contigs) num_threads(nthreads)
	for (size_t i = 0; i < cont_num; i++ ) {

		//if ( ordered_contigs[i].first > 20000)
		INFO("processing contig" << ordered_contigs[i].second << " of length " << ordered_contigs[i].first <<" in thread number " << omp_get_thread_num());
//		ContigProcessor pc(all_contigs_copy[ordered_contigs[i].second].sam_filenames[0].first, all_contigs_copy[ordered_contigs[i].second].input_contig_filename);
		ContigProcessor pc(all_contigs_copy[ordered_contigs[i].second].sam_filenames, all_contigs_copy[ordered_contigs[i].second].input_contig_filename);

		pc.process_multiple_sam_files();
		//pc.process_all_sam_files();
//		INFO("contig " << ordered_contigs[i].second << " in thread number " << omp_get_thread_num() <<" processed");
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
