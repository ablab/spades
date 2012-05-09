#pragma once

namespace compare {
template<size_t k>
void FillBagForStrand(const Sequence& strand,
		map<Seq<k>, size_t, typename Seq<k>::less2>& bag) {
	if (strand.size() < k)
		return;
	Seq<k> kmer(strand);
	kmer >> 'A';
	for (size_t i = k - 1; i < strand.size(); ++i) {
		kmer = kmer << strand[i];
		bag[kmer] += 1;
	}
}

template<size_t k>
void FillRepeats(const Sequence& genome,
		set<Seq<k>, typename Seq<k>::less2>& repeats) {
	map<Seq<k>, size_t, typename Seq<k>::less2> bag;

	FillBagForStrand(genome, bag);
	FillBagForStrand(!genome, bag);

	for (auto it = bag.begin(); it != bag.end(); ++it) {
		if (it->second > 1)
			repeats.insert(it->first);
	}
}

//todo
struct ClearedGenome {
	typedef vector<MappingRange> Ranges;
	Sequence s;
	Ranges mapping;
};

//template<size_t k>
//class RepeatCleaner {
//	typedef Seq<k> Kmer;
//	typedef set<Kmer, typename Kmer::less2> Repeats;
//
//	void MarkPositions(size_t start, size_t end, vector<bool>& answer) {
//		for (size_t i = start; i < end; ++i) {
//			answer[i] = true;
//		}
//	}
//
//	void MarkRepeatNucls(const Sequence& s, const Repeats& repeats, vector<bool>& answer) {
////		vector<bool> answer(s.size(), false);
//		Kmer kmer(s);
//		kmer = kmer >> 'A';
//		for (size_t i = k - 1 ; i < s.size(); ++i) {
//			kmer = kmer << s[i];
//			if (repeats.count(kmer) > 0) {
//				MarkPositions(i - k + 1, i + 1, answer);
//			}
//		}
//	}
//
//	void MarkShortIslands(, size_t threshold = k) {
//
//	}
//
//public:
//
//};

template<size_t k>
Sequence ClearGenome(const Sequence& genome,
		const set<Seq<k>, typename Seq<k>::less2>& repeats) {
	INFO("Clearing genome");
	if (genome.size() < k)
		return genome;

	string answer;
	for (size_t i = 0; i < k - 1; ++i) {
		answer += nucl(genome[i]);
	}
	//intervals of kmers that should be resolved afterwards
	vector<Range> repeat_intervals;
	Seq<k> kmer(genome);
	size_t curr_pos = 0;
	//curr_pos + k - next nucl pos
	bool changed = false;
	while (curr_pos + k != genome.size()) {
		size_t int_start = curr_pos;
		while (repeats.count(kmer) > 0 && curr_pos + k < genome.size()) {
			kmer = kmer << genome[curr_pos + k];
			curr_pos++;
			changed = true;
		}
		if (int_start != curr_pos)
			repeat_intervals.push_back(Range(int_start, curr_pos));

		if (curr_pos + k == genome.size())
			break;

		while (repeats.count(kmer) == 0 && curr_pos + k < genome.size()) {
			answer += nucl(kmer[k - 1]);
			kmer = kmer << genome[curr_pos + k];
			curr_pos++;
		}
	}
	if (changed) {
		INFO("Genome was changed during cleaning");
	} else {
		INFO("Genome wasn't changed during cleaning");
	}
	return Sequence(answer);
}

template<size_t k>
Sequence ClearGenome(const Sequence& genome) {
	INFO("Clearing genome of repeats");

	set<Seq<k>, typename Seq<k>::less2> repeats;
	INFO("Filling set of repeats");
	FillRepeats<k>(genome, repeats);
	INFO("Clearing genome");
	return ClearGenome<k>(genome, repeats);
}

//todo bad strategy for assembly cleaning
template<size_t k>
pair<Sequence, vector<Sequence>> Clear(const Sequence& genome,
		const vector<Sequence>& assembly) {
	INFO("Clearing genome of repeats");

	set<Seq<k>, typename Seq<k>::less2> repeats;
	INFO("Filling set of repeats");
	FillRepeats<k>(genome, repeats);
	for (auto it = assembly.begin(); it != assembly.end(); ++it) {
		FillRepeats(*it, repeats);
	}INFO("Clearing genome");
	Sequence new_genome = ClearGenome<k>(genome, repeats);
	INFO("Clearing assembly");
	vector<Sequence> new_assembly;
	for (auto it = assembly.begin(); it != assembly.end(); ++it) {
		new_assembly.push_back(ClearGenome<k>(*it, repeats));
	}
	return make_pair(new_genome, new_assembly);
}

template<size_t k>
pair<Sequence, Sequence> ClearGenomes(const pair<Sequence, Sequence>& genomes) {
	INFO("Clearing genomes from repeats");

	set<Seq<k>, typename Seq<k>::less2> repeats;
	INFO("Filling set of repeats");
	FillRepeats<k>(genomes.first, repeats);
	FillRepeats<k>(genomes.second, repeats);
	INFO("Clearing genomes");
	return make_pair(ClearGenome<k>(genomes.first, repeats),
			ClearGenome<k>(genomes.second, repeats));
}

template<size_t k>
pair<Sequence, Sequence> TotallyClearGenomes(
		const pair<Sequence, Sequence>& genomes) {
	static const size_t iter_count = 1;
	pair<Sequence, Sequence> tmp = genomes;
	for (size_t i = 0; i < iter_count; ++i) {
		INFO("Cleaning iteration " << i);
		tmp = ClearGenomes<k>(tmp);
	}
	return tmp;
}

template<size_t k>
bool CheckNoRepeats(const Sequence& genome) {
	set<Seq<k>, typename Seq<k>::less2> repeats;
	FillRepeats<k>(genome, repeats);
	return repeats.empty();
}

}
