//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once


#include "standard.hpp"
#include "logger/logger.hpp"
#include "de/paired_info.hpp"
#include "omni/omni_utils.hpp"
#include "omni/abstract_conjugate_graph.hpp"
#include "omni/abstract_nonconjugate_graph.hpp"
#include "utils.hpp"
#include "debruijn_graph.hpp"
#include "levenshtein.hpp"

#include "omni/omni_tools.hpp"

#include "omni/id_track_handler.hpp"
#include "omni/edge_labels_handler.hpp"
#include "omni/graph_component.hpp"

class MySamRecord{
public:
	string QNAME;
	size_t FLAG;
	string RNAME;
	size_t POS;
	size_t MAPQ;
	string CIGAR;
	string RNEXT;
	size_t PNEXT;
	int TLEN;
	string SEQ;
	string QUAL;
	MySamRecord(string QNAME_){
		//Default constructor produce unallign record
		QNAME = QNAME_;
		FLAG = 4;
		RNAME = "*";
		POS = 0;
		MAPQ = 255;
		CIGAR = "*";
		RNEXT = "*";
		PNEXT = 0;
		TLEN = 0;
		SEQ = "*";
		QUAL = "*";
	}

	MySamRecord(string QNAME_, string SEQ_, string QUAL_ = "*"){
		//Default constructor produce unalign record
		QNAME = QNAME_;
		FLAG = 4;
		RNAME = "*";
		POS = 0;
		MAPQ = 255;
		CIGAR = "*";
		RNEXT = "*";
		PNEXT = 0;
		TLEN = 0;
		SEQ = SEQ_;
		QUAL = QUAL_;
//		for (size_t i = 0; i<SEQ.size(); i++) QUAL = QUAL+"I";
	}
	string str(){
		return QNAME+"\t"+ToString(FLAG)+"\t"+RNAME + "\t"+ ToString(POS)+ "\t"+ToString(MAPQ) +"\t"
				+ CIGAR +"\t"+	RNEXT + "\t" + ToString(PNEXT) +"\t" + ToString(TLEN) + "\t"+ SEQ + "\t"+ QUAL;

	}
	string map_str(){
		if (FLAG & 0x10)
			return QNAME+"\t-\t" + RNAME + "\t" + ToString(POS) + "\t"+ SEQ + "\t"+ QUAL + "\t0";
		else
			return QNAME+"\t+\t" + RNAME + "\t" + ToString(POS) + "\t"+ SEQ + "\t"+ QUAL + "\t0";
	}
	bool is_aligned(){
		return !(FLAG & 0x04);
	}
	bool is_rc(){
		return (FLAG & 0x10);
	}
	void updateInfoFromNextSegment(MySamRecord & NextSegment){
		PNEXT = NextSegment.POS;
		if (RNAME == NextSegment.RNAME){
			RNEXT = "=";
		}
		else RNEXT = NextSegment.RNAME;
		if (NextSegment.FLAG & 0x10)
			FLAG |= 0x20;
	}
};


template<class ComparingObject>
int CountDiference(ComparingObject const& s_r, ComparingObject const& orig_s_r, int shift){
	if (s_r.size() + shift > orig_s_r.size()){
		return orig_s_r.size();
	}
	int diff = 0;
	for (int i = 0; i < (int)s_r.size(); i++){
		if (s_r[i] != orig_s_r[i + shift] ) diff++;
	}
	return diff;
}

size_t CountDiference(const io::SingleRead& s_r, const io::SingleRead& orig_s_r, size_t shift){
	if (s_r.size() + shift > orig_s_r.size()){
		return orig_s_r.size();
	}
	size_t diff = 0;
	for (size_t i = 0; i < s_r.size(); i++){
		if (s_r[i] != orig_s_r[i + shift] ) diff++;
	}
	return diff;
}
pair<size_t, size_t> SubstitutionShifts(const io::SingleRead& s_r, const io::SingleRead& orig_s_r, double threshhold ){
	size_t difference = orig_s_r.size() - s_r.size();
	size_t best_i = 0;
	size_t best_diff = orig_s_r.size();
	for (size_t i=0; i <= difference; i++){
		size_t cur_diff = CountDiference(s_r, orig_s_r, i);
		if (math::gr((double)(s_r.size() - cur_diff), threshhold * (double)s_r.size())) return (make_pair(i,difference - i));
		else {
			if (cur_diff < best_diff) {
				best_diff = cur_diff;
				best_i = i;
			}
		}
	}
	return (make_pair(best_i,difference - best_i));
}
void SubstituteByOriginalRead(MySamRecord& MySam, const io::SingleRead& s_r, const io::SingleRead& orig_s_r, bool print_quality = false){
	io::SingleRead orig = orig_s_r;
	if (MySam.FLAG & 0x10) orig = !orig;
	if (s_r.size() < orig.size()){
		pair<size_t, size_t> shifts = SubstitutionShifts(s_r, orig_s_r, 0.9);
//		pair<int, int> shifts = s_r.position_in_original();
//		shifts.second = orig_s_r.size() - shifts.second;
		MySam.SEQ = orig.sequence().str();
		if (print_quality)
			MySam.QUAL = orig.GetPhredQualityString(); //orig.GetPhredQualityString()
		else MySam.QUAL = "*";

		if (MySam.FLAG & 0x10) MySam.CIGAR = (shifts.second != 0 ? ToString(shifts.second) + "S": "") + MySam.CIGAR + (shifts.first != 0 ? ToString(shifts.first) + "S": "");
		else 				   MySam.CIGAR = (shifts.first != 0 ? ToString(shifts.first) + "S": "") + MySam.CIGAR + (shifts.second != 0 ? ToString(shifts.second) + "S": "");
	}
	else
	{
		MySam.SEQ = orig.sequence().str();
		if (print_quality)
			MySam.QUAL = orig.GetPhredQualityString();
		else MySam.QUAL = "*";
	}
}






template<class Graph, class SequenceMapper>
class BaseInternalAligner {
protected:

    size_t k_;

	typedef typename Graph::EdgeId EdgeId;
	const Graph &graph_;
	const SequenceMapper& mapper_;
	bool adjust;
	bool map_mode;
	bool print_broken;
	bool print_quality;
	map<EdgeId, pair<string, bool>> SeqNames;
	FILE* samOut;
	size_t ProcessedReads;
	size_t SplittedReads;
	size_t SuccesfullReads;
	size_t SamRecordsCount;

	string  ProduceNodeName(EdgeId edge, size_t id){
		size_t len = graph_.length(edge) + k_; //to do: check it
		double cov = graph_.coverage(edge);
		return "NODE_" + ToString(id) + "_length_" +ToString(len)  + "_cov_" + ToString(cov);
	}

	void FillRefSequences() {
		size_t id = 0;
		for (auto it = graph_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			if (SeqNames.find(*it) == SeqNames.end()) {
				string s = ProduceNodeName(*it, id);
				id++;
				if (! map_mode) fprintf(samOut, "@SQ\tSN:%s\tLN:%u\n", s.c_str(), (unsigned int)(graph_.length(*it) + k_));
				SeqNames.insert(make_pair(*it, make_pair(s,true)));
				EdgeId conj = graph_.conjugate(*it);
				if (conj != *it)
					SeqNames.insert(make_pair(conj, make_pair(s,false)));
			}
		}
		DEBUG("Seq names prepared");
	}

	void InitializeSamFile(const string& sam_output_filename){
		INFO("STAGE == SAM file printing");
		samOut = fopen(sam_output_filename.c_str(), "w");
 		if (samOut) {
			if (! map_mode) fprintf(samOut, "@HD\tVN:1.4\tSO:unsorted\n");
			FillRefSequences();
		}
		ProcessedReads = 0;
		SuccesfullReads = 0;
		SplittedReads = 0;
		SamRecordsCount = 0;

	}
	void FinalizeSamFile(){
		if (samOut) {
			fclose(samOut);
			samOut = NULL;
		}
		DEBUG("SAM file generating finished");
		INFO("Processed "<< ProcessedReads<< " reads, "<< SuccesfullReads<< " ("<< ((double)SuccesfullReads*100)/ProcessedReads<<"%) of them aligned");
//		INFO( SplittedReads<<  " ( "<< ((double)SplittedReads*100)/ProcessedReads<<"%)"<< " reads split on few edges.");
//		INFO( " SamRecordsCount "<< SamRecordsCount);
	}

	MySamRecord CreateSingleSAMFromRange(string read_name, Sequence read, EdgeId edge, MappingRange range, bool rc, string qual = "*"){
		Sequence ref_seq = this->graph_.EdgeNucls(edge);
		size_t ref_start = range.mapped_range.start_pos;
		size_t ref_end = range.mapped_range.end_pos+k_;

		pair<pair<int, int>, string> cigar_pair;
		if (this->adjust) {
			if (CountDiference<Sequence>(read.Subseq(range.initial_range.start_pos, range.initial_range.end_pos + k_), ref_seq.Subseq(ref_start, ref_end), 0) != 0)
			{
				if (ref_start > 10) ref_start -= 10;
				else ref_start = 0;

				if (ref_end + 10 < ref_seq.size()) ref_end += 10;
				else ref_end = ref_seq.size();

				cigar_pair = best_edit_distance_cigar(read.Subseq(range.initial_range.start_pos, range.initial_range.end_pos+k_).str()
																				, ref_seq.Subseq(ref_start, ref_end).str());
				ref_start = cigar_pair.first.first + ref_start+1;
			}
			else {
				cigar_pair = make_pair(make_pair(0, 0), ToString((int)(range.initial_range.end_pos+k_ - range.initial_range.start_pos))+"M");
				ref_start = range.mapped_range.start_pos+1;

			}
		} else
		{
			cigar_pair = make_pair(make_pair(0, 0), ToString((int)(range.initial_range.end_pos+k_ - range.initial_range.start_pos))+"M");
			ref_start = range.mapped_range.start_pos + 1;
		}
//		MySamRecord SamRec(read_name, read.Subseq(range.initial_range.start_pos, range.initial_range.end_pos+k_).str());
		MySamRecord SamRec(read_name, read.str(), qual);
		SamRec.FLAG = 0;
		SamRec.POS = ref_start; //path1[0].second.mapped_range.start_pos+1;
		SamRec.RNAME = this->SeqNames[edge].first;
		SamRec.CIGAR = "";
		if (range.initial_range.start_pos > 0)
			SamRec.CIGAR = ToString((int)(range.initial_range.start_pos))+"S";
		SamRec.CIGAR += cigar_pair.second;
		if (read.size() > range.initial_range.end_pos+k_)
			SamRec.CIGAR += ToString((int)(read.size() - (range.initial_range.end_pos+k_)))+"S";
		if (rc) SamRec.FLAG |= 0x10;
//			if (path1[0].second.mapped_range.end_pos - path1[0].second.mapped_range.start_pos != path1[0].second.initial_range.end_pos - path1[0].second.initial_range.start_pos){
//				WARN("Possible indels: "<< s_r.name()<<" VS "<< SamRec.RNAME);
//			}
		return SamRec;

	}

	void UpdateAlignedPair(MySamRecord& SamRec1, MySamRecord& SamRec2){
		size_t found;
		found = SamRec1.QNAME.rfind("/");
		if (found!=string::npos)
			SamRec1.QNAME.replace (found,2,"");
		SamRec2.QNAME = SamRec1.QNAME;
		SamRec1.FLAG |= 0x01;// | 0x40;
		SamRec2.FLAG |= 0x01;// | 0x80;
		SamRec1.FLAG |= 0x02;
		SamRec2.FLAG |= 0x02;
		SamRec1.updateInfoFromNextSegment(SamRec2);
		SamRec2.updateInfoFromNextSegment(SamRec1);
		if ((SamRec1.RNEXT == "=")&&(SamRec1.POS > SamRec1.PNEXT)){
			MySamRecord TmpSam = SamRec1;
			SamRec1 = SamRec2;
			SamRec2 = TmpSam;
		}
		if (SamRec1.RNEXT == "=") {
			SamRec1.TLEN = (int)SamRec1.PNEXT + (int)SamRec2.SEQ.size() - (int)SamRec1.POS;
			SamRec2.TLEN = -SamRec1.TLEN;
		}
		SamRec1.FLAG |= 0x40;
		SamRec2.FLAG |= 0x80;
	}
	void UpdatePair(MySamRecord& SamRec1, MySamRecord& SamRec2){
		if (SamRec1.is_aligned() && SamRec2.is_aligned()) {
			UpdateAlignedPair(SamRec1, SamRec2);
		}
		else {
			size_t found;
			found = SamRec1.QNAME.rfind("/");
			if (found!=string::npos)
				SamRec1.QNAME.replace (found,2,"");
			SamRec2.QNAME = SamRec1.QNAME;
			SamRec1.FLAG |= 0x01;// | 0x40;
			SamRec2.FLAG |= 0x01;// | 0x80;
			SamRec1.updateInfoFromNextSegment(SamRec2);
			SamRec2.updateInfoFromNextSegment(SamRec1);
			if (!SamRec1.is_aligned()) SamRec2.FLAG |= 0x08;
			if (!SamRec2.is_aligned()) SamRec1.FLAG |= 0x08;
			SamRec1.FLAG |= 0x40;
			SamRec2.FLAG |= 0x80;
		}
	}



	virtual void ProcessSingleRead(const io::SingleRead& s_r) {
//		if (map_mode){
//			MySamRecord SamRec = CreateSingleSAMFromSingleRead(s_r);
//			if (SamRec.is_aligned())
//				fprintf(samOut, "%s\n", SamRec.map_str().c_str());
//		}
//		else {
//			fprintf(samOut, "%s\n", CreateSingleSAMFromSingleRead(s_r).str().c_str());
//		}
	}

	virtual void ProcessPairedRead(const io::PairedRead& p_r) {
	}


public:
	BaseInternalAligner(size_t k, Graph& g, SequenceMapper& m, bool adjust_reads = false, bool output_map_format = false, bool print_broken_pairs = false, bool print_quality_ = false): k_(k),
		graph_(g), mapper_(m), adjust(adjust_reads), map_mode(output_map_format), print_broken(print_broken_pairs), print_quality(print_quality_)
	{};

	template<class Stream>
	void AlignSingleReads(Stream& s, const string& sam_output_filename) {
		InitializeSamFile(sam_output_filename);
		while (!s.eof()) {
			io::SingleRead s_r;
			s >> s_r;
			ProcessSingleRead(s_r);
		}
		FinalizeSamFile();


	}
	template<class Stream>
	void AlignPairedReads(Stream& s, const string& sam_output_filename) {
		size_t n = 0;
		InitializeSamFile(sam_output_filename);
		while (!s.eof()) {
			io::PairedRead p_r;
			s >> p_r;
			VERBOSE_POWER(++n, " paired reads processed");
			ProcessPairedRead(p_r);
		}
		FinalizeSamFile();
	}

	template<class PairedStream, class SingleStream>
	void AlignReads(PairedStream& p_s, SingleStream& s_s, const string& sam_output_filename) {
		InitializeSamFile(sam_output_filename);
		while (!p_s.eof()) {
				io::PairedRead p_r;
				p_s >> p_r;
				ProcessPairedRead(p_r);
		}
		while (!s_s.eof()) {
				io::SingleRead s_r;
				s_s >> s_r;
				ProcessSingleRead(s_r);
		}
		FinalizeSamFile();
	}



};



template<class Graph, class SequenceMapper>
class SimpleInternalAligner: public BaseInternalAligner<Graph, SequenceMapper>{
protected:
	typedef typename Graph::EdgeId EdgeId;

	MySamRecord CreateSingleSAMFromSingleRead(const io::SingleRead& s_r, bool try_rc_align_first = false){
		Sequence read = s_r.sequence();
		if (try_rc_align_first) read = !read;
		MappingPath<EdgeId> path1 = this->mapper_.MapSequence(read);
		this->ProcessedReads++;
		string qual = "*";
		if (path1.size() == 1){
			EdgeId edge = path1[0].first;
			bool rc = try_rc_align_first;
			if (!this->SeqNames[edge].second){
				read = !read;
				path1 = this->mapper_.MapSequence(read);
				edge = path1[0].first;
				rc = !try_rc_align_first;
				if (this->print_quality) {
					qual = (!s_r).GetPhredQualityString();
				}
			}
			else {
				if (this->print_quality) {
					qual = s_r.GetPhredQualityString();
				}
			}
			this->SuccesfullReads++;
			return this->CreateSingleSAMFromRange(s_r.original_name(), read, edge, path1[0].second, rc, qual);

		} else {
			if (this->print_quality) {
				qual = s_r.GetPhredQualityString();
			}
			if (path1.size() > 1) this->SplittedReads++;
			MySamRecord SamRec(s_r.original_name(), s_r.GetSequenceString(), qual);
			return SamRec;
		}
	}


	void ProcessSingleRead(const io::SingleRead& s_r) {
		if (this->map_mode){
			MySamRecord SamRec = CreateSingleSAMFromSingleRead(s_r);
			if (SamRec.is_aligned())
				fprintf(this->samOut, "%s\n", SamRec.map_str().c_str());
		}
		else {
			fprintf(this->samOut, "%s\n", CreateSingleSAMFromSingleRead(s_r).str().c_str());
		}
	}

	void ProcessPairedRead(const io::PairedRead& p_r) {
		MySamRecord SamRec1 = CreateSingleSAMFromSingleRead(p_r.first());
		MySamRecord SamRec2 = CreateSingleSAMFromSingleRead(p_r.second(), !SamRec1.is_rc());
		if (this->map_mode) {
			if (SamRec1.is_aligned())
				fprintf(this->samOut, "%s\n", SamRec1.map_str().c_str());
			if (SamRec2.is_aligned())
				fprintf(this->samOut, "%s\n", SamRec2.map_str().c_str());
		}
		else {
			if (SamRec1.is_aligned()&&SamRec2.is_aligned()){
				if (!this->print_broken) {
					this->UpdateAlignedPair(SamRec1, SamRec2);
					fprintf(this->samOut, "%s\n", SamRec1.str().c_str());
					fprintf(this->samOut, "%s\n", SamRec2.str().c_str());
				}
			}
			if (this->print_broken) {
				this->UpdatePair(SamRec1, SamRec2);
				fprintf(this->samOut, "%s\n", SamRec1.str().c_str());
				fprintf(this->samOut, "%s\n", SamRec2.str().c_str());
			}
		}
	}

	void PrepareHeader();

public:
	SimpleInternalAligner(size_t k, Graph& g, SequenceMapper& m, bool adjust_reads = false, bool output_map_format = false, bool print_broken_pairs = false, bool print_quality = false):
		BaseInternalAligner<Graph, SequenceMapper>(k, g, m, adjust_reads, output_map_format, print_broken_pairs, print_quality)
//		graph_(g), mapper_(m), adjust(adjust_reads), map_mode(output_map_format), print_broken(print_broken_pairs)
	{};
};



template<class Graph, class SequenceMapper>
class ResolvedInternalAligner: public BaseInternalAligner<Graph, SequenceMapper>{
protected:
	typedef typename Graph::EdgeId EdgeId;
	const Graph &original_graph_;
	EdgeLabelHandler<Graph> &convertor_;


//	vector<MySamRecord> CreateMultipleSAMFromSingleRead(const io::SingleRead& s_r){
//	}
	inline Range ReverceRanges(Range& range, size_t length){
//		INFO("Reverce range " << range<< " with len "<< length);
		return Range(length - range.end_pos, length - range.start_pos);
	}

	Range ConvertRanges(const Range& range, EdgeId& edge, vector<EdgeId> labels, bool rc) {
		size_t shift = 0;
		size_t i = 0;
		while ((i < labels.size() && labels[i] != edge) ){
			shift += original_graph_.length(labels[i]);
			i++;
		}

//		INFO("Range "<<range<<" shifted on "<<shift);
		if (i < labels.size()) {
			if (rc) {
				size_t whole_length = shift;
				while ((i < labels.size()) ){
					whole_length += original_graph_.length(labels[i]);
					i++;
				}
				DEBUG("Range is rc, whole length" << whole_length);
				return Range(whole_length - range.end_pos - shift, whole_length - range.start_pos - shift);
			} else {
				return Range(range.start_pos + shift, range.end_pos + shift);
			}
		}
		else {
			ERROR("Can not convert range");
			return Range(0,0);
		}
	}

	Range ConvertRanges(MappingPath<EdgeId>& path, vector<EdgeId> labels, bool rc) {
		size_t start_shift = 0;
		size_t end_shift = 0;
		size_t converted = 0;
		size_t i = 0;
		size_t start_i = 0;

//		for(size_t j = 0; j<path.size();  j++){
//			INFO(" " << path[j].first<<" "<<path[j].second.mapped_range);
//		}
//		INFO("labs "<<labels);

		while (i < labels.size()) {
			if ((labels[i] != path[0].first) ){
				start_shift += original_graph_.length(labels[i]);
				end_shift += original_graph_.length(labels[i]);
				i++;
				start_i++;
			}
			else {
				break;
			}
		}

		if (i < labels.size()) {
//			converted++;
//			i++;
			while ((converted<path.size()) && (i < labels.size())){
//				INFO("i "<<i<<" conv "<<converted);
				if (labels[i] == path[converted].first) {
					converted++;
					if (converted<path.size()) {
						end_shift += original_graph_.length(labels[i]);
						i++;
					}
				} else {
					converted = 0;
					start_shift += original_graph_.length(labels[start_i]);
					end_shift = start_shift;
					i = start_i + 1;
					start_i++;
//					INFO("bad i "<<i<<" conv "<<converted);
					while (i < labels.size()) {
						if ((labels[i] != path[0].first) ){
							start_shift += original_graph_.length(labels[i]);
							end_shift += original_graph_.length(labels[i]);
							i++;
							start_i++;
						} else {
							break;
						}
					}
				}
			}
		}

		if (i < labels.size()){
//			INFO("done i "<<i<<" conv "<<converted);
			if (rc) {
				size_t whole_length = end_shift;
				while ((i < labels.size()) ){
					whole_length += original_graph_.length(labels[i]);
					i++;
				}
//				INFO("Range is rc, whole length" << whole_length);
				return Range(whole_length - path[path.size()-1].second.mapped_range.end_pos - end_shift, whole_length - path[0].second.mapped_range.start_pos - start_shift);
			} else {
//				INFO("Start shift " << start_shift<<" end shift "<<end_shift);
				return Range(path[0].second.mapped_range.start_pos + start_shift, path[path.size()-1].second.mapped_range.end_pos + end_shift);
			}
		}
		else {

//			ERROR("Can not convert range");
//			for(size_t j = 0; j<path.size();  j++){
//				INFO(" " << path[j].first<<" "<<path[j].second.mapped_range);
//			}
//			INFO("labs "<<labels);
			return Range(0,0);
		}
	}

	vector<MySamRecord> CreateMultipleSAMFromSingleRead(const io::SingleRead& s_r, bool try_rc_align_first = false){
		vector<MySamRecord> result;
		Sequence proto_read = s_r.sequence();
		if (try_rc_align_first) proto_read = !proto_read;
		MappingPath<EdgeId> path1 = this->mapper_.MapSequence(proto_read);
		this->ProcessedReads++;
		if (path1.size() > 0){
			EdgeId proto_edge = path1[0].first;

//			set<VertexId> my_set;
			for (auto iter = convertor_.edge_inclusions[proto_edge].begin(); iter != convertor_.edge_inclusions[proto_edge].end(); ++iter){
				bool rc = false;
				EdgeId edge = *iter;
				Range i_r(path1[0].second.initial_range.start_pos, path1[path1.size()-1].second.initial_range.end_pos);
				Sequence read;
				string qual = "*";
				if (!this->SeqNames[edge].second){
					read = !proto_read;
					rc = true;
					i_r = ReverceRanges(i_r, read.size() - BaseInternalAligner<Graph, SequenceMapper>::k_);
					if (this->print_quality) {
						qual = (!s_r).GetPhredQualityString();
					}

				} else {
					read = proto_read;
					if (this->print_quality) {
						qual = s_r.GetPhredQualityString();
					}
				}
//				Range m_r = ConvertRanges(path1[0].second.mapped_range, proto_edge, convertor_.edge_labels[edge], rc);
				Range m_r = ConvertRanges(path1, convertor_.edge_labels[edge], rc);
				MappingRange new_range(i_r, m_r);
				if (rc) edge = this->graph_.conjugate(edge);
				if (m_r.end_pos!=0) {
					result.push_back(this->CreateSingleSAMFromRange(s_r.original_name(), read, edge, new_range, try_rc_align_first?!rc:rc, qual));
					this->SamRecordsCount++;
				}
			}
			if (result.size()==0) {
				result.push_back(MySamRecord(s_r.original_name(), s_r.GetSequenceString(), (this->print_quality ? s_r.GetPhredQualityString():"*")));
				this->SamRecordsCount++;
			}
			else {
				this->SuccesfullReads++;
			}


		} else {
//			if (path1.size() > 1) {
//				this->SplittedReads++;
//			}
			MySamRecord SamRec(s_r.original_name(), s_r.GetSequenceString(), (this->print_quality ? s_r.GetPhredQualityString():"*"));
			this->SamRecordsCount++;
			result.push_back(SamRec);
		}
		return result;
	}


public:
	ResolvedInternalAligner(size_t k, Graph& g, Graph& orig_g, SequenceMapper& m, EdgeLabelHandler<Graph>& EdgeConversionHandler, bool adjust_reads = false, bool output_map_format = false, bool print_broken_pairs = false, bool print_quality = false):
		BaseInternalAligner<Graph, SequenceMapper>(k, g, m, adjust_reads, output_map_format, print_broken_pairs, print_quality), original_graph_(orig_g), convertor_(EdgeConversionHandler)
//		graph_(g), mapper_(m), adjust(adjust_reads), map_mode(output_map_format), print_broken(print_broken_pairs)
	{};

	void ProcessSingleRead(const io::SingleRead& s_r) {
		vector<MySamRecord> SamRecs = CreateMultipleSAMFromSingleRead(s_r);
		if (this->map_mode){
			for(size_t i = 0; i < SamRecs.size(); i++){
				if (SamRecs[i].is_aligned())
					fprintf(this->samOut, "%s\n", SamRecs[i].map_str().c_str());
			}
		}
		else {
			for(size_t i = 0; i < SamRecs.size(); i++){
				fprintf(this->samOut, "%s\n", SamRecs[i].str().c_str());
			}
		}
	}


	void ProcessPairedRead(const io::PairedRead& p_r) {
		if (this->map_mode) {
			ProcessSingleRead(p_r.first());
			ProcessSingleRead(p_r.second());
		}
		else {
			vector<MySamRecord> SamRecs1 = CreateMultipleSAMFromSingleRead(p_r.first());
			bool try_rc_first = false;
			if (SamRecs1.size() > 0){
				if (!SamRecs1[0].is_rc()) try_rc_first = true;
			}
			vector<MySamRecord> SamRecs2 = CreateMultipleSAMFromSingleRead(p_r.second(), try_rc_first);
//			INFO( " SamRecordsCount "<< this->SamRecordsCount);
			if ((SamRecs1.size() == 1)&&((SamRecs1.size() == 1))){
				if (SamRecs1[0].is_aligned()&&SamRecs2[0].is_aligned()){
					if (!this->print_broken) {
						this->UpdateAlignedPair(SamRecs1[0], SamRecs2[0]);
						fprintf(this->samOut, "%s\n", SamRecs1[0].str().c_str());
						fprintf(this->samOut, "%s\n", SamRecs2[0].str().c_str());
					}
				}
				if (this->print_broken) {
					this->UpdatePair(SamRecs1[0], SamRecs2[0]);
					fprintf(this->samOut, "%s\n", SamRecs1[0].str().c_str());
					fprintf(this->samOut, "%s\n", SamRecs2[0].str().c_str());
				}
			}
		}
	}
};




template<class Graph, class SequenceMapper>
class OriginalReadsResolvedInternalAligner : public ResolvedInternalAligner<Graph, SequenceMapper>{
protected:
//	void SubstituteByOriginalRead(MySamRecord& MySam, const io::SingleRead& s_r, const io::SingleRead& orig_s_r ){
//		io::SingleRead orig = orig_s_r;
//		if (MySam.FLAG & 0x10) orig = !orig;
//		if (s_r.size() < orig.size()){
//			MySam.SEQ = orig.sequence().str();
//			MySam.QUAL = orig.GetPhredQualityString();
//			if (MySam.FLAG & 0x10) MySam.CIGAR = ToString((int)(orig_s_r.size() - s_r.size())) + "S"+MySam.CIGAR;
//			else MySam.CIGAR += ToString((int)(orig_s_r.size() - s_r.size())) + "S";
//		}
//		else
//		{
//			MySam.SEQ = orig.sequence().str();
//			MySam.QUAL = orig.GetPhredQualityString();
//		}
//	}

public:
	OriginalReadsResolvedInternalAligner(size_t k, Graph& g, Graph& orig_g, SequenceMapper& m, EdgeLabelHandler<Graph>& EdgeConversionHandler, bool adjust_reads = false, bool output_map_format = false, bool print_broken_pairs = false, bool print_quality = false):
		ResolvedInternalAligner<Graph, SequenceMapper>(k, g, orig_g, m, EdgeConversionHandler, adjust_reads, output_map_format, print_broken_pairs, print_quality)
	{};

	template<class Stream, class OrigStream>
	void AlignPairedReads(OrigStream& original_s, Stream& s, const string& sam_output_filename) {
		this->InitializeSamFile(sam_output_filename);
		size_t n = 0;
		while (!s.eof()){
			io::PairedRead p_r;
			s >> p_r;
			while (!original_s.eof()) {
				io::PairedRead orig_p_r;
				original_s >> orig_p_r;
				if (p_r.first().original_name() == orig_p_r.first().original_name().substr(0, p_r.first().original_name().size())){
					VERBOSE_POWER(++n, " paired reads processed");
					ProcessPairedReadWithOriginal(p_r, orig_p_r);
					break;
				}
			}
		}
		this->FinalizeSamFile();
	}

	void ProcessPairedReadWithOriginal(const io::PairedRead& p_r, const io::PairedRead& orig_p_r ) {
		if (p_r.first().size()>orig_p_r.first().size()) return;
		if (p_r.second().size()>orig_p_r.second().size()) return;
		if (this->map_mode) {
			ProcessSingleReadWithOriginal(p_r.first(), orig_p_r.first());
			ProcessSingleReadWithOriginal(p_r.second(), orig_p_r.second());
		}
		else {
			vector<MySamRecord> SamRecs1 = this->CreateMultipleSAMFromSingleRead(p_r.first());
			bool try_rc_first = false;
			if (SamRecs1.size() > 0){
				if (!SamRecs1[0].is_rc()) try_rc_first = true;
			}
			vector<MySamRecord> SamRecs2 = this->CreateMultipleSAMFromSingleRead(p_r.second(), try_rc_first);
//			INFO( " SamRecordsCount "<< this->SamRecordsCount);
			if ((SamRecs1.size() == 1)&&((SamRecs1.size() == 1))){
				SubstituteByOriginalRead(SamRecs1[0], p_r.first(), orig_p_r.first(), this->print_quality);
				SubstituteByOriginalRead(SamRecs2[0], p_r.second(), orig_p_r.second(), this->print_quality);
				if (SamRecs1[0].is_aligned()&&SamRecs2[0].is_aligned()){
					if (!this->print_broken) {
						this->UpdateAlignedPair(SamRecs1[0], SamRecs2[0]);
						fprintf(this->samOut, "%s\n", SamRecs1[0].str().c_str());
						fprintf(this->samOut, "%s\n", SamRecs2[0].str().c_str());
					}
				}
				if (this->print_broken) {
					this->UpdatePair(SamRecs1[0], SamRecs2[0]);
					fprintf(this->samOut, "%s\n", SamRecs1[0].str().c_str());
					fprintf(this->samOut, "%s\n", SamRecs2[0].str().c_str());
				}
			}
		}
	}


	void ProcessSingleReadWithOriginal(const io::SingleRead& s_r, const io::SingleRead& orig_s_r) {
		if (s_r.size()>orig_s_r.size()) return;
		vector<MySamRecord> SamRecs = this->CreateMultipleSAMFromSingleRead(s_r);
		if (this->map_mode){
			for(size_t i = 0; i < SamRecs.size(); i++){
				if (SamRecs[i].is_aligned()){
					SubstituteByOriginalRead(SamRecs[i], s_r, orig_s_r, this->print_quality);
					fprintf(this->samOut, "%s\n", SamRecs[i].map_str().c_str());
				}
			}
		}
		else {
			for(size_t i = 0; i < SamRecs.size(); i++){
				SubstituteByOriginalRead(SamRecs[i], s_r, orig_s_r, this->print_quality);
				fprintf(this->samOut, "%s\n", SamRecs[i].str().c_str());
			}
		}
	}
};


template<class Graph, class SequenceMapper>
class OriginalReadsSimpleInternalAligner : public SimpleInternalAligner<Graph, SequenceMapper>{
protected:
public:
	OriginalReadsSimpleInternalAligner(size_t k, Graph& g, SequenceMapper& m, bool adjust_reads = false, bool output_map_format = false, bool print_broken_pairs = false, bool print_quality = false):
		SimpleInternalAligner<Graph, SequenceMapper>(k, g, m, adjust_reads, output_map_format, print_broken_pairs, print_quality)
	{};

	template<class Stream, class OrigStream>
	void AlignPairedReads(OrigStream& original_s, Stream& s, const string& sam_output_filename) {
		this->InitializeSamFile(sam_output_filename);
		size_t n = 0;
		while (!s.eof()){
			io::PairedRead p_r;
			s >> p_r;
			while (!original_s.eof()) {
				io::PairedRead orig_p_r;
				original_s >> orig_p_r;
				if (p_r.first().original_name() == orig_p_r.first().original_name().substr(0, p_r.first().original_name().size())){
					VERBOSE_POWER(++n, " paired reads processed");
					ProcessPairedReadWithOriginal(p_r, orig_p_r);
					break;
				}
			}
		}
		this->FinalizeSamFile();
	}

	template<class Stream, class OrigStream>
	void AlignSingleReads(OrigStream& original_s, Stream& s, const string& sam_output_filename) {
		this->InitializeSamFile(sam_output_filename);
		size_t n = 0;
		while (!s.eof()){
			io::SingleRead s_r;
			s >> s_r;
			while (!original_s.eof()) {
				io::SingleRead orig_s_r;
				original_s >> orig_s_r;
				if (s_r.original_name() == orig_s_r.original_name().substr(0, s_r.original_name().size())){
					VERBOSE_POWER(++n, " single reads processed");
					ProcessSingleReadWithOriginal(s_r, orig_s_r);
					break;
				}
			}
		}
		this->FinalizeSamFile();
	}


	void ProcessPairedReadWithOriginal(const io::PairedRead& p_r, const io::PairedRead& orig_p_r ) {
		if (p_r.first().size()>orig_p_r.first().size()) return;
		if (p_r.second().size()>orig_p_r.second().size()) return;

		MySamRecord SamRec1 = this->CreateSingleSAMFromSingleRead(p_r.first());
		MySamRecord SamRec2 = this->CreateSingleSAMFromSingleRead(p_r.second(), !SamRec1.is_rc());
		if (this->map_mode) {
			SubstituteByOriginalRead(SamRec1, p_r.first(), orig_p_r.first(), this->print_quality);
			SubstituteByOriginalRead(SamRec2, p_r.second(), orig_p_r.second(), this->print_quality);
			if (SamRec1.is_aligned())
				fprintf(this->samOut, "%s\n", SamRec1.map_str().c_str());
			if (SamRec2.is_aligned())
				fprintf(this->samOut, "%s\n", SamRec2.map_str().c_str());
		}
		else {
			SubstituteByOriginalRead(SamRec1, p_r.first(), orig_p_r.first(), this->print_quality);
			SubstituteByOriginalRead(SamRec2, p_r.second(), orig_p_r.second(), this->print_quality);

			if (SamRec1.is_aligned()&&SamRec2.is_aligned()){
				if (!this->print_broken) {
					this->UpdateAlignedPair(SamRec1, SamRec2);
					fprintf(this->samOut, "%s\n", SamRec1.str().c_str());
					fprintf(this->samOut, "%s\n", SamRec2.str().c_str());
				}
			}
			if (this->print_broken) {
				this->UpdatePair(SamRec1, SamRec2);
				fprintf(this->samOut, "%s\n", SamRec1.str().c_str());
				fprintf(this->samOut, "%s\n", SamRec2.str().c_str());
			}
		}
	}





	void ProcessSingleReadWithOriginal(const io::SingleRead& s_r, const io::SingleRead& orig_s_r) {
		if (s_r.size()>orig_s_r.size()) return;
		if (this->map_mode){
			MySamRecord SamRec = this->CreateSingleSAMFromSingleRead(s_r);
			if (SamRec.is_aligned()){
				SubstituteByOriginalRead(SamRec, s_r, orig_s_r, this->print_quality);
				fprintf(this->samOut, "%s\n", SamRec.map_str().c_str());
			}
		}
		else {
			MySamRecord SamRec = this->CreateSingleSAMFromSingleRead(s_r);
			SubstituteByOriginalRead(SamRec, s_r, orig_s_r, this->print_quality);
			fprintf(this->samOut, "%s\n", SamRec.str().c_str());
		}
	}



};


