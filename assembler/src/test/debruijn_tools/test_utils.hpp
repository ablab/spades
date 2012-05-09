#pragma once

namespace compare {

inline double uniform_01() {
	static boost::mt19937 rng(43);
	static boost::uniform_01<boost::mt19937> zeroone(rng);
	return zeroone();
}

inline bool event_happened(double rate) {
	return ls(uniform_01(), rate);
}

inline int rand_int(size_t min, size_t max) {
	static boost::mt19937 rng(43);
	boost::uniform_int<> un_int(min, max);
	boost::variate_generator<boost::mt19937&, boost::uniform_int<> > die(rng,
			un_int);
	return die();
}

inline char switch_nucl(char n) {
	VERIFY(is_nucl(n));
	return nucl((dignucl(n) + rand_int(1, 3)) % 4);
}

inline Sequence IntroduceReversal(const Sequence& s, size_t min_len, size_t max_len) {
	VERIFY(s.size() > min_len);
	//inclusive
	size_t start = rand_int(0, s.size() - min_len);
	size_t len = rand_int(min_len, std::min(max_len, s.size() - start));
	//exclusive
	size_t end = start + len;
	INFO(
			"Reversing fragment of length " << len << " from " << start << " to " << end);
	return s.Subseq(0, start) + !s.Subseq(start, end) + s.Subseq(end);
}

inline Sequence IntroduceReversals(const Sequence& s, size_t rev_count, size_t min_len,
		size_t max_len) {
	Sequence res = s;
	for (size_t i = 0; i < rev_count; ++i) {
		res = IntroduceReversal(res, min_len, max_len);
	}
	return res;
}

inline Sequence IntroduceMutations(const Sequence& s, double rate) {
	VERIFY(ge(rate, 0.) && ls(rate, 1.0));
	string as_str = s.str();
	for (size_t i = 0; i < s.size(); ++i) {
		if (event_happened(rate)) {
			as_str[i] = switch_nucl(as_str[i]);
		}
	}
	return Sequence(as_str);
}

template<class gp_t>
void ConstructRepeatGraph(gp_t& gp) {
	io::VectorReader<io::SingleRead> stream(
			io::SingleRead("genome", gp.genome.str()));
	io::RCReaderWrapper<io::SingleRead> rc_stream(stream);
	ConstructGraph<gp_t::k_value, typename gp_t::graph_t>(gp.g, gp.index,
			rc_stream);
}

template<class Graph>
vector<Sequence> EdgesSequences(const Graph& g) {
	vector<Sequence> res;
	set<EdgeId> edges;
	for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
		if (edges.find(*it) == edges.end()) {
			res.push_back(g.EdgeNucls(*it));
			edges.insert(g.conjugate(*it));
		}
	}
	return res;
}

template<class gp_t>
vector<Sequence> RepeatGraphEdges(const Sequence& genome) {
	typedef typename gp_t::graph_t Graph;
	typedef typename Graph::EdgeId EdgeId;

	gp_t gp(genome);
	ConstructRepeatGraph(gp);
	return EdgesSequences(gp.g);
}


bool CheckFileDiff(const string& file1, const string& file2) {
	INFO("Checking differences between " << file1 << " and " << file2);
	checkFileExistenceFATAL(file1);
	checkFileExistenceFATAL(file2);
	ifstream f1(file1.c_str());
	ifstream f2(file2.c_str());
	while (!f1.eof() && !f2.eof()) {
		string s1;
		f1 >> s1;
		string s2;
		f2 >> s2;
		if (s1 != s2)
			return false;
	}
	if (!f1.eof() || !f2.eof())
		return false;
	return true;
}

}
