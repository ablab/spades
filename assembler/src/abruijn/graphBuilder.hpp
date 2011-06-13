#ifndef GRAPHBUILDER_H_
#define GRAPHBUILDER_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <ext/hash_map>
#include "hash.hpp"
#include "parameters.hpp"
#include "logging.hpp"
#include "omnigraph.hpp"
#include "ireadstream.hpp"

namespace abruijn {

using namespace std;
using namespace __gnu_cxx;
using hashing::hash_t;

typedef omnigraph::Omnigraph Graph;

class GraphBuilder
{
	//typedef hash_map< Sequence, int, HashSym<Sequence>, EqSym<Sequence> > SeqCount;
	//SeqCount seqCount;
	/**
	 * Stores bitmask:
	 * 1 = the lexicographically lesser k-mer has been seen as non-rightmost k-mer in read
	 * 2 = the greater k-mer ...
	 */
	hashing::HashSym<Sequence> hashSym;
	typedef vector<hash_t> hash_vector;
	hash_vector ha;

public:
	set<hash_t> earmarked_hashes;
	typedef map<hash_t, char> RightExtensions;
	RightExtensions has_right;
	RightExtensions tips;
	typedef map<hash_t, size_t> TipExtenstionTable;
	typedef map<hash_t, TipExtenstionTable> TipExtensions;
	TipExtensions tip_extensions;
	hash_vector hbest;
//	typedef abruijn::Graph Graph;
	Graph graph_;
//	typedef hash_map <Sequence, Graph::VertexId, hashing::HashSym<Sequence>, hashing::EqSym<Sequence> > SeqVertice;
	typedef hash_map <Sequence, Graph::VertexId, hashing::Hash<Sequence> > SeqVertice;
	SeqVertice seqVertice;
	size_t htake_;

	void findMinimizers(Sequence s);
	void findLocalMinimizers(Sequence s, size_t window_size);
	void findSecondMinimizer(Sequence s);
	void revealTips(Sequence s);
	void findTipExtensions(Sequence s);
	void lookRight(Sequence s);
	void addToGraph(Sequence s);
	bool hasVertex(Sequence s);
	Graph::VertexId getOrCreateVertex(Sequence s);
	Graph::VertexId createVertex(Sequence s);

private:
	DECL_LOGGER("GraphBuilder")
};

template <typename Reader>
class GraphBuilderMaster {
	Reader reader_;
	int mode_;
	GraphBuilder gb_;
public:
	GraphBuilderMaster(Reader reader, size_t htake, int mode) :
		reader_(reader), mode_(mode), gb_() {
		gb_.hbest.reserve(htake);
		gb_.htake_ = htake;
	}

	void build() {
		Read r;

		INFO("===== Finding " << gb_.htake_ << " minimizers in each read... =====");
		reader_.reset();
		for (size_t i = 0; !reader_.eof(); ++i) {
			reader_ >> r;
			if (mode_ & 1) {
				gb_.findMinimizers(r.getSequence());
			} else {
				gb_.findLocalMinimizers(r.getSequence(), 51);
			}
			VERBOSE(i, " single reads");
		}
		INFO("Done: " << gb_.earmarked_hashes.size() << " earmarked hashes");

		if ((mode_ & 2) && (gb_.htake_ == 1)) {
			INFO("===== Finding second minimizers... =====");
			reader_.reset();
			for (size_t i = 0; !reader_.eof(); ++i) {
				reader_ >> r;
				gb_.findSecondMinimizer(r.getSequence());
				VERBOSE(i, " single reads");
			}
			INFO("Done: " << gb_.earmarked_hashes.size() << " earmarked hashes");
		}

		for (int tip_iteration = 0;; tip_iteration++) {
			size_t eh = gb_.earmarked_hashes.size();
			gb_.has_right.clear();
			gb_.tips.clear();
			gb_.tip_extensions.clear();

			INFO("===== Revealing tips... =====");
			reader_.reset();
			for (size_t i = 0; !reader_.eof(); ++i) {
				reader_ >> r;
				gb_.revealTips(r.getSequence());
				VERBOSE(i, " single reads");
			}
			for (map<hash_t, char>::iterator it = gb_.has_right.begin(); it != gb_.has_right.end(); ++it) {
				if (it->second != 3) {
					TRACE(it->first << " " << (int) it->second);
					gb_.tips.insert(*it);
				}
			}
			gb_.has_right.clear();
			INFO("Done: " << gb_.tips.size() << " tips.");
			if (gb_.tips.size() == 0) {
				INFO("Quitting tip extension procedure");
				break;
			}

			INFO("===== Finding tip extensions... =====");
			reader_.reset();
			for (size_t i = 0; !reader_.eof(); ++i) {
				reader_ >> r;
				gb_.findTipExtensions(r.getSequence());
				VERBOSE(i, " single reads");
			}
			INFO("Done: " << gb_.has_right.size() << " possible tip extensions");
			if (gb_.has_right.size() == 0) {
				INFO("Quitting tip extension procedure");
				break;
			}

			INFO("===== Looking to the right... =====");
			reader_.reset();
			for (size_t i = 0; !reader_.eof(); ++i) {
				reader_ >> r;
				gb_.lookRight(r.getSequence());
				VERBOSE(i, " single reads");
			}
			for (GraphBuilder::TipExtensions::iterator it = gb_.tip_extensions.begin(); it != gb_.tip_extensions.end(); ++it) {
				TRACE("Trying to extend tip " << it->first);
				bool ok = false;
				hash_t found = hashing::kMax;
				size_t best_dist = 0;
				for (GraphBuilder::TipExtenstionTable::iterator ext = it->second.begin(); ext != it->second.end(); ext++) {
					if (gb_.earmarked_hashes.count(ext->first)) {
						ok = true;
						TRACE("Wonderful hit!");
						break;
					}
					if (gb_.has_right[ext->first] == 3) {
						TRACE("Candidate at distance " << ext->second);
						found = ext->first;
						best_dist = -1;
						continue;
					}
					if (gb_.has_right[ext->first] > 0) {
						TRACE("Dubious candidate at distance " << ext->second);
						if (ext->second > best_dist) {
							found = ext->first;
							best_dist = ext->second;
						}
						continue;
					}
				}
				if (ok) {
					continue;
				}
				TRACE("Selected candidate: " << best_dist);
				gb_.earmarked_hashes.insert(found);
			}
			INFO("Done: " << eh << " -> " << gb_.earmarked_hashes.size() << " earmarked hashes");
			if (eh == gb_.earmarked_hashes.size()) {
				INFO("Quitting tip extension procedure");
				break;
			}
		}

		INFO("===== Adding reads to graph as paths... =====");
		reader_.reset();
		for (size_t i = 0; !reader_.eof(); ++i) {
			reader_ >> r;
			gb_.addToGraph(r.getSequence());
			VERBOSE(i, " single reads");
		}
		INFO("Done: " << gb_.graph_.size() << " vertices");

//		INFO("===== Condensing-A graph... =====");
//		gb_.graph_.Condense();
//		INFO("Done: " << gb_.graph_.vertices.size() << " vertices");

		return;
	}

	omnigraph::Omnigraph* graph() {
		return &gb_.graph_;
	}

	void SpellGenomeThroughGraph (size_t cut) {
		/// we assume here that the graph is already built

		/// outputting A Bruijn graph
//		ofstream outfile("./data/abruijn/abruijnspellgenome.dot", ios::out);
//		gvis::GraphPrinter<Vertex*> printer("abruijnspellgenome", outfile);
//		for (Vertices::iterator v = graph()->vertices.begin(); v != graph()->vertices.end(); ++v) {
//			printer.addVertex(*v, ToString(**v));
//			for (Edges::iterator it = (*v)->edges().begin(); it != (*v)->edges().end(); ++it) {
//				printer.addEdge(*v, it->first, ToString(it->second));
//			}
//		}


		/// reading the reference genome
		string const ref_genome_filename = "./data/input/MG1655-K12.fasta.gz";
		Read ref_genome;

		ireadstream genome_stream(ref_genome_filename);
		genome_stream >> ref_genome;
		genome_stream.close();

		Sequence ref_seq = ref_genome.getSequence();
		if ( cut > 0 )
			ref_seq = ref_seq.Subseq(0, cut);

		//INFO ( "ref_seq: " << ref_seq << " " << ref_seq.size() );

		/// computing hash-values of all the K-mers of the reference genome
		hashing::HashSym<Sequence> hashsym;
		vector<hash_t> ha;
		ha.resize(ref_seq.size()-K+1);
		hashsym.kmers(ref_seq, ha);

		size_t num_of_missing_kmers = 0;
		size_t num_of_earmarked_kmers = 0;
		size_t num_of_missing_edges = 0;
		size_t num_of_missing_lengths = 0;

		int previous_index = -1, current_index = -1;
		Sequence previous_kmer ( "" ), current_kmer ( "" );

		for (unsigned int i = 0; i != ha.size (); ++i ) {
			if (gb_.earmarked_hashes.count(ha[i])) {
				INFO(i << " out of " << ha.size () << " (" << i*100 / ha.size () << "%)" );
				++num_of_earmarked_kmers;

				current_index = i;
				current_kmer  = ref_genome.getSequence().Subseq(i,i+K);
				if (!gb_.hasVertex( current_kmer)) {
					++num_of_missing_kmers;
					INFO("k-mer " << current_kmer << " is present in the genome, but not in the graph");
					continue;
				}

				if (-1 == previous_index) {
					previous_index = i;
					previous_kmer  = current_kmer;
					//printer.threadStart( graph()->getVertex( previous_kmer ), 4*graph()->vertices.size() );
				}
				else {
					current_index = i;
					current_kmer  = ref_genome.getSequence().Subseq( i,i+K );
					//printer.threadAdd( graph()->getVertex( current_kmer) );

					size_t edge_length = current_index - previous_index;
					assert ( edge_length > 0 );

					assert (gb_.hasVertex(previous_kmer) && gb_.hasVertex(current_kmer));
					Graph::VertexId previous_vertex = gb_.getOrCreateVertex(previous_kmer);
					Graph::VertexId current_vertex  = gb_.getOrCreateVertex(current_kmer);
					assert ( previous_vertex && current_vertex );

					vector<Graph::EdgeId> edges = graph()->GetEdgesBetween(previous_vertex, current_vertex);
					if (edges.size() == 0) {
						INFO ( "missing edge from " << previous_kmer << " to " << current_kmer );
						++num_of_missing_edges;
					}
					else {
						bool occ = false;
						for (auto it = edges.begin(); it != edges.end(); ++it) {
							occ |= (graph()->data(*it).length() == edge_length);
						}
						if (!occ) {
							INFO ( "missing length" );
							++num_of_missing_lengths;
						}
					}

					previous_index = current_index;
					previous_kmer  = current_kmer;
				}
			} // if earmarked
		} // for

//		printer.output();
//		outfile.close ();

		/// printing stats
		INFO ( "number of earmarked k-mers: " << num_of_earmarked_kmers );
		INFO ( "number of missing k-mers: " << num_of_missing_kmers );
		INFO ( "number of missing edges: " << num_of_missing_edges );
		INFO ( "number of missing lengths: " << num_of_missing_lengths );

//		return ( num_of_missing_kmers + num_of_missing_edges + num_of_missing_lengths == 0 );
		return;
	}

private:
	DECL_LOGGER("GraphBuildMaster")
};

} // namespace abruijn

#endif /* GRAPHBUILDER_H_ */
