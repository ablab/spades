#pragma once


#include "standard.hpp"
#include "logging.hpp"
#include "omni/paired_info.hpp"
#include "omni/omni_utils.hpp"
#include "omni/abstract_conjugate_graph.hpp"
#include "omni/abstract_nonconjugate_graph.hpp"
#include "utils.hpp"
#include "new_debruijn.hpp"
#include "levenshtein.hpp"

#include "omni/omni_tools.hpp"
#include "omni/omnigraph.hpp"

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

	MySamRecord(string QNAME_, string SEQ_){
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
		SEQ = SEQ_;
//		QUAL = "*";
		for (size_t i = 0; i<SEQ.size(); i++) QUAL = QUAL+"I";
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


template<class Graph, class SequenceMapper>
class BaseInternalAligner {
protected:
	typedef typename Graph::EdgeId EdgeId;
	const Graph &graph_;
	const SequenceMapper& mapper_;
	bool adjust;
	bool map_mode;
	bool print_broken;
	map<EdgeId, pair<string, bool>> SeqNames;
	FILE* samOut;
	size_t ProcessedReads;
	size_t SplittedReads;
	size_t SuccesfullReads;
	size_t SamRecordsCount;

	string  ProduceNodeName(EdgeId edge, size_t id){
		size_t len = graph_.length(edge) + debruijn_graph::K; //to do: check it
		double cov = graph_.coverage(edge);
		return "NODE_" + ToString(id) + "_length_" +ToString(len)  + "_cov_" + ToString(cov);
	}

	void FillRefSequences() {
		size_t id = 0;
		for (auto it = graph_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			if (SeqNames.find(*it) == SeqNames.end()) {
				string s = ProduceNodeName(*it, id);
				id++;
				if (! map_mode) fprintf(samOut, "@SQ\tSN:%s\tLN:%u\n", s.c_str(), (unsigned int)(graph_.length(*it) + debruijn_graph::K));
				SeqNames.insert(make_pair(*it, make_pair(s,true)));
				EdgeId conj = graph_.conjugate(*it);
				if (conj != *it)
					SeqNames.insert(make_pair(conj, make_pair(s,false)));
			}
		}INFO("Seq names prepared");
	}

	void InitializeSamFile(const string& sam_output_filename){
		INFO("SAM file initialize");
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
		INFO("SAM file generating finished");
		INFO("Processed "<< ProcessedReads<< " reads, "<< SuccesfullReads<< " ( "<< ((double)SuccesfullReads*100)/ProcessedReads<<"%) of them aligned.");
		INFO( SplittedReads<<  " ( "<< ((double)SplittedReads*100)/ProcessedReads<<"%)"<< " reads split on few edges.");
		INFO( " SamRecordsCount "<< SamRecordsCount);
	}

	MySamRecord CreateSingleSAMFromRange(string read_name, Sequence read, EdgeId edge, MappingRange range, bool rc){
		Sequence ref_seq = this->graph_.EdgeNucls(edge);
		size_t ref_start = range.mapped_range.start_pos;
		size_t ref_end = range.mapped_range.end_pos+debruijn_graph::K;
		if (ref_start > 10) ref_start -= 10;
		else ref_start = 0;

		if (ref_end + 10 < ref_seq.size()) ref_end += 10;
		else ref_end = ref_seq.size();

		pair<pair<int, int>, string> cigar_pair;
		if (this->adjust) {
//			INFO(" init "<<range.initial_range<<" map "<< range.mapped_range);
			cigar_pair = best_edit_distance_cigar(read.Subseq(range.initial_range.start_pos, range.initial_range.end_pos+debruijn_graph::K).str()
																			, ref_seq.Subseq(ref_start, ref_end).str());
//			INFO(" init "<<cigar_pair);
			ref_start = cigar_pair.first.first + ref_start+1;
		} else
		{
			cigar_pair = make_pair(make_pair(0, 0), ToString((int)(range.initial_range.end_pos+debruijn_graph::K - range.initial_range.start_pos))+"M");
			ref_start = range.mapped_range.start_pos;
		}
		MySamRecord SamRec(read_name, read.Subseq(range.initial_range.start_pos, range.initial_range.end_pos+debruijn_graph::K).str());
		SamRec.FLAG = 0;
		this->SuccesfullReads++;
		SamRec.POS = ref_start; //path1[0].second.mapped_range.start_pos+1;
		SamRec.RNAME = this->SeqNames[edge].first;
		SamRec.CIGAR = cigar_pair.second;
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
	BaseInternalAligner(Graph& g, SequenceMapper& m, bool adjust_reads = false, bool output_map_format = false, bool print_broken_pairs = false):
		graph_(g), mapper_(m), adjust(adjust_reads), map_mode(output_map_format), print_broken(print_broken_pairs)
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
		InitializeSamFile(sam_output_filename);
		while (!s.eof()) {
			io::PairedRead p_r;
			s >> p_r;
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
private:
	typedef typename Graph::EdgeId EdgeId;

	MySamRecord CreateSingleSAMFromSingleRead(const io::SingleRead& s_r){
		Sequence read = s_r.sequence();
		MappingPath<EdgeId> path1 = this->mapper_.MapSequence(read);
		this->ProcessedReads++;
		if (path1.size() == 1){
			EdgeId edge = path1[0].first;
			bool rc = false;
			if (!this->SeqNames[edge].second){
				read = !read;
				path1 = this->mapper_.MapSequence(read);
				edge = path1[0].first;
				rc = true;
			}
			return CreateSingleSAMFromRange(s_r.name(), read, edge, path1[0].second, rc);

		} else {
			if (path1.size() > 1) this->SplittedReads++;
			MySamRecord SamRec(s_r.name(), s_r.GetSequenceString());
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
		MySamRecord SamRec2 = CreateSingleSAMFromSingleRead(p_r.second());
		if (this->map_mode) {
			if (SamRec1.is_aligned())
				fprintf(this->samOut, "%s\n", SamRec1.map_str().c_str());
			if (SamRec2.is_aligned())
				fprintf(this->samOut, "%s\n", SamRec2.map_str().c_str());
		}
		else {
			if (SamRec1.is_aligned()&&SamRec2.is_aligned()){
				this->UpdateAlignedPair(SamRec1, SamRec2);
				if (!this->print_broken) {
					fprintf(this->samOut, "%s\n", SamRec1.str().c_str());
					fprintf(this->samOut, "%s\n", SamRec2.str().c_str());
				}
			}
			if (this->print_broken) {
				fprintf(this->samOut, "%s\n", SamRec1.str().c_str());
				fprintf(this->samOut, "%s\n", SamRec2.str().c_str());
			}
		}
	}

	void PrepareHeader();

public:
	SimpleInternalAligner(Graph& g, SequenceMapper& m, bool adjust_reads = false, bool output_map_format = false, bool print_broken_pairs = false):
		BaseInternalAligner<Graph, SequenceMapper>(g, m, adjust_reads, output_map_format, print_broken_pairs)
//		graph_(g), mapper_(m), adjust(adjust_reads), map_mode(output_map_format), print_broken(print_broken_pairs)
	{};
};



template<class Graph, class SequenceMapper>
class ResolvedInternalAligner: public BaseInternalAligner<Graph, SequenceMapper>{
private:
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

			ERROR("Can not convert range");
					for(size_t j = 0; j<path.size();  j++){
						INFO(" " << path[j].first<<" "<<path[j].second.mapped_range);
					}
					INFO("labs "<<labels);
			return Range(0,0);
		}
	}

	vector<MySamRecord> CreateMultipleSAMFromSingleRead(const io::SingleRead& s_r){
		vector<MySamRecord> result;
		Sequence proto_read = s_r.sequence();
		MappingPath<EdgeId> path1 = this->mapper_.MapSequence(proto_read);
		this->ProcessedReads++;
		if (path1.size() > 0){
			EdgeId proto_edge = path1[0].first;
			for (auto iter = convertor_.edge_inclusions[proto_edge].begin(); iter != convertor_.edge_inclusions[proto_edge].end(); ++iter){
				bool rc = false;
				EdgeId edge = *iter;
				Range i_r(path1[0].second.initial_range.start_pos, path1[path1.size()-1].second.initial_range.end_pos);
				Sequence read;
				if (!this->SeqNames[edge].second){
					read = !proto_read;
					rc = true;
					i_r = ReverceRanges(i_r, read.size() - debruijn_graph::K);
				} else {
					read = proto_read;
				}
//				Range m_r = ConvertRanges(path1[0].second.mapped_range, proto_edge, convertor_.edge_labels[edge], rc);
				Range m_r = ConvertRanges(path1, convertor_.edge_labels[edge], rc);
				MappingRange new_range(i_r, m_r);
				if (rc) edge = this->graph_.conjugate(edge);
				if (m_r.end_pos!=0) {
					result.push_back(CreateSingleSAMFromRange(s_r.name(), read, edge, new_range, rc));
					this->SamRecordsCount++;
				}
			}
			if (result.size()==0) {
				result.push_back(MySamRecord(s_r.name(), s_r.GetSequenceString()));
				this->SamRecordsCount++;
			}


		} else {
//			if (path1.size() > 1) {
//				this->SplittedReads++;
//			}
			MySamRecord SamRec(s_r.name(), s_r.GetSequenceString());
			this->SamRecordsCount++;
			result.push_back(SamRec);
		}
		return result;
	}


public:
	ResolvedInternalAligner(Graph& g, Graph& orig_g, SequenceMapper& m, EdgeLabelHandler<Graph>& EdgeConversionHandler, bool adjust_reads = false, bool output_map_format = false, bool print_broken_pairs = false):
		BaseInternalAligner<Graph, SequenceMapper>(g, m, adjust_reads, output_map_format, print_broken_pairs), original_graph_(orig_g), convertor_(EdgeConversionHandler)
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
			vector<MySamRecord> SamRecs2 = CreateMultipleSAMFromSingleRead(p_r.second());
//			INFO( " SamRecordsCount "<< this->SamRecordsCount);
			if ((SamRecs1.size() == 1)&&((SamRecs1.size() == 1))){
				if (SamRecs1[0].is_aligned()&&SamRecs2[0].is_aligned()){
					this->UpdateAlignedPair(SamRecs1[0], SamRecs2[0]);
					if (!this->print_broken) {
						fprintf(this->samOut, "%s\n", SamRecs1[0].str().c_str());
						fprintf(this->samOut, "%s\n", SamRecs2[0].str().c_str());
					}
				}
				if (this->print_broken) {
					fprintf(this->samOut, "%s\n", SamRecs1[0].str().c_str());
					fprintf(this->samOut, "%s\n", SamRecs2[0].str().c_str());
				}
			}
		}
	}
};

