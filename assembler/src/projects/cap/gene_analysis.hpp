//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/standard_base.hpp"
#include "utils/simple_tools.hpp"
#include "comparison_utils.hpp"
#include "boost/tokenizer.hpp"
#include "coloring.hpp"

//todo deprecated
namespace cap {
using namespace omnigraph;

typedef Sequence Genome;
typedef map<Genome, Range> GeneCoordinates;

//range of nucleotide positions on main strand; true if main strand
typedef pair<Range, bool> Pos;
typedef size_t GenomeId;
typedef size_t GeneId;

typedef multimap<GenomeId, Pos> GenePositions;
typedef vector<Range> Coordinates;

struct GenomeInfo {
  GenomeId id;
  string name;
  string short_name;
  Sequence sequence;

  GenomeInfo()
      : id(0) {

  }

  GenomeInfo(const GenomeId& id_, const string& name_,
             const string& short_name_, const Sequence& sequence_)
      : id(id_),
        name(name_),
        short_name(short_name_),
        sequence(sequence_) {
  }

};

struct GeneInfo {
  GeneId id;
//    string name;
  GenePositions gene_positions;

  GeneInfo()
      : id(0) {
  }

  GeneInfo(const GeneId& id_/*, const string& name_*/)
      : id(id_)/*, name(name_)*/{
  }

  void AddPosition(const GenomeId& genome_id, const Pos& pos) {
    gene_positions.insert(make_pair(genome_id, pos));
  }

};

Range VerifiedRange(Range r, size_t k, size_t genome_length) {
  VERIFY(genome_length > k + 1);
  VERIFY(r.start_pos < genome_length);
  VERIFY(r.end_pos <= genome_length);
  VERIFY(r.start_pos <= r.end_pos);
  return r;
}

//ALL k-mers, lying in the nucl range
Range NuclToKmerRange(Range nucl_range, size_t k, size_t genome_length) {
//    return Range(std::max(0, int(nucl_range.start_pos) - int(k)),
//                 nucl_range.end_pos);

  size_t start_pos =
      (nucl_range.start_pos + k + 1 > genome_length) ?
          genome_length - k - 1 : nucl_range.start_pos;
  size_t end_pos = std::max((int) nucl_range.end_pos - (int) k,
                            int(start_pos + 1));

  return VerifiedRange(Range(start_pos, end_pos), k, genome_length);
}

//ALL nucls, covered by k-mer range
Range KmerToNuclRange(Range kmer_range, size_t k, size_t genome_length) {
//    return Range(min(kmer_range.start_pos + k, kmer_range.end_pos - 1),
//                 kmer_range.end_pos);
  return VerifiedRange(Range(kmer_range.start_pos, kmer_range.end_pos + k + 1),
                       k, genome_length);
}

//todo fix
Range OppositeStrandNuclCoord(Range nucl_coord, size_t genome_length) {
  VERIFY(nucl_coord.end_pos <= genome_length);
  return Range(genome_length - 1 - nucl_coord.end_pos,
               genome_length - 1 - nucl_coord.start_pos);
}

//Updates k-mer coordinates to the coordinates of start/end of condensed edge path in simplified graph
// (gene is somewhere inside it)
template<class gp_t>
class CoordinatesUpdater {
  typedef typename gp_t::graph_t Graph;
  typedef typename Graph::EdgeId EdgeId;
  typedef typename Graph::VertexId VertexId;

  const gp_t& gp_;

  //Works with and returns k-mer coordinates!
  Range NewCoords(const MappingPath<EdgeId>& path, Range coord) const {
    size_t cumm_length = 0;
    Range answer(0, 0);
    size_t i = 0;
    for (;
        i < path.size()
            && path[i].second.initial_range.end_pos <= coord.start_pos; ++i) {
      cumm_length += gp_.g.length(path[i].first);
    }
    VERIFY(i < path.size());
    VERIFY(path[i].second.initial_range.end_pos > coord.start_pos);
    VERIFY(path[i].second.initial_range.start_pos <= coord.start_pos);
    answer.start_pos = cumm_length;
    for (;
        i < path.size() && path[i].second.initial_range.end_pos < coord.end_pos;
        ++i) {
      cumm_length += gp_.g.length(path[i].first);
    }
    VERIFY(i < path.size());
    VERIFY(path[i].second.initial_range.end_pos >= coord.end_pos);
    VERIFY(path[i].second.initial_range.start_pos < coord.end_pos);
    answer.end_pos = cumm_length + gp_.g.length(path[i].first);
    return answer;
  }

 public:
  CoordinatesUpdater(const gp_t& gp)
      : gp_(gp) {

  }

  pair<Genome, Coordinates> Update(const Genome& genome,
                                   const Coordinates& coords) const {
    pair<Genome, Coordinates> answer;
    auto mapper = *MapperInstance(gp_);
    auto mapping_path = mapper.MapSequence(genome);
    for (Range r : coords) {
      answer.second.push_back(
          KmerToNuclRange(
              NewCoords(
                  mapping_path,
                  NuclToKmerRange(r, gp_.k_value,
                      genome.size())),
              gp_.k_value, genome.size()));
    }
    auto read_corrector = GraphReadCorrectorInstance(gp_.g, mapper);
    answer.first = read_corrector->Modify(genome);
    return answer;
  }
 private:
  DECL_LOGGER("CoordinatesUpdater")
  ;
};

//deals with nucleotide coordinates!
struct GeneCollection {
  map<GenomeId, GenomeInfo> genomes;
  map<GeneId, GeneInfo> genes;

  map<string, GenomeId> genome_name_id_mapping;
 private:

  void LoadGenomes(const set<string>& genome_names,
                   const string& genomes_folder) {
    size_t id = 0;
    for (string name : genome_names) {
    string filename = genomes_folder + name;
    path::CheckFileExistenceFATAL(filename);
    genomes.insert(
        make_pair(
            id,
            GenomeInfo(id, name, name, ReadGenome(filename))));
    genome_name_id_mapping.insert(make_pair(name, id));
    id++;
  }
}

GenomeId genome_id(const string& name) const {
  return get(genome_name_id_mapping, name);
}

void LoadGenomes(const string& file_with_genomes,
    const string& genomes_folder) {
  path::CheckFileExistenceFATAL(file_with_genomes);
  ifstream stream(file_with_genomes);
  set<string> genome_names;
  string name;
  while (!stream.eof()) {
    stream >> name;
    genome_names.insert(name);
  }
  LoadGenomes(genome_names, genomes_folder);
}

void SaveGenomes(const string& /* folder */) const {
  for (auto it = genomes.begin(); it != genomes.end(); ++it) {
    WriteGenome(it->second.sequence, it->second.name);
  }
}

void SaveGeneInfo(const string& filename) const {
  ofstream stream(filename);
  for (auto it1 = genes.begin(); it1 != genes.end(); ++it1) {
    for (auto it = it1->second.gene_positions.begin(); it != it1->second.gene_positions.end(); ++it) {
      stream << (boost::format("%i\t%s\t%i\t%i\t%s\n") % it1->second.id
          % get(genomes, it->first).name
          % it->second.first.start_pos
          % it->second.first.end_pos
          % it->second.second).str();
    }
  }
}

set<int> LoadGeneIDs(const string& file_with_ids) {
  path::CheckFileExistenceFATAL(file_with_ids);
  ifstream stream(file_with_ids);
  set<int> gene_ids;
  int id;
  while (!stream.eof()) {
    stream >> id;
    gene_ids.insert(id);
  }
  return gene_ids;
}

void AddGeneInfo(const GeneId& gene_id, const GenomeId& genome_id, const Range& range, bool forward) {
  if (genes.count(gene_id) == 0) {
    genes.insert(make_pair(gene_id, GeneInfo(gene_id)));
  }
  get(genes, gene_id).AddPosition(genome_id, Pos(range, forward));
}

//ortholog ids is better wording
void LoadGeneInfo(const string& filename, set<int> gene_ids) {
  using boost::tokenizer;
  using boost::escaped_list_separator;
  path::CheckFileExistenceFATAL(filename);
  ifstream stream(filename);
  string line;
  while (!stream.eof()) {
    stream >> line;
    tokenizer<escaped_list_separator<char>> tk(
        line, escaped_list_separator<char>('\t'));
    vector<string> record(tk.begin(), tk.end());
    //0 - id
    //1 - genome name
    //2 - start
    //3 - end
    //4 - forward/reverse
    VERIFY(record.size() == 8);
    VERIFY(record[4] == "reverse" || record[4] == "forward");
    int gene_id = lexical_cast<int>(record[0]);
    if (gene_ids.count(gene_id) > 0) {
      AddGeneInfo(gene_id, genome_id(record[1]),
          Range(lexical_cast<size_t>(record[2]), lexical_cast<size_t>(record[3]))
          , record[4] == "forward");
    }
  }
}

public:

GeneCollection() {}

//ortholog ids is better wording
void Load(const string& file_with_genomes, const string& genomes_folder,
    const string& file_with_gene_info, const string& file_with_ids) {
  LoadGenomes(file_with_genomes, genomes_folder);
  LoadGeneInfo(file_with_gene_info, LoadGeneIDs(file_with_ids));
}

void Save(const string& root_folder, const string& genomes_folder,
    const string& file_with_gene_info) const {
  SaveGenomes(root_folder + genomes_folder);
  SaveGeneInfo(root_folder + file_with_gene_info);
}

template<class gp_t>
void Update(const gp_t& gp) {
  CoordinatesUpdater<gp_t> updater(gp);
  for (GenomeId genome_id : key_set(genomes)) {
    Coordinates gene_coords;
    vector<GeneId> gene_ids;
    for (GeneId gene_id : key_set(genes)) {
      for (Pos pos:  get_all(genes[gene_id].gene_positions, genome_id)) {
        gene_ids.push_back(gene_id);
        VERIFY(pos.second);
        gene_coords.push_back(pos.first);
      }
    }
    auto updated = updater.Update(genomes.find(genome_id)->second.sequence, gene_coords);
    genomes.find(genome_id)->second.sequence = updated.first;

    //clearing gene positions
    for (GeneId gene_id : key_set(genes)) {
      genes[gene_id].gene_positions.clear();
    }

    //updating gene positions
    for (size_t j = 0; j < gene_ids.size(); ++j) {
      genes[gene_ids[j]].gene_positions.insert(
          make_pair(genome_id,
              make_pair(updated.second[j], true)));
    }
  }
}

private:
DECL_LOGGER("GeneCollection")
;
};

//template<class gp_t>
//void WriteGeneLocality(const GeneCollection& gene_collection, const gp_t& gp,
//                       const string& folder,
//                       const ColorHandler<typename gp_t::graph_t>& coloring) {
//  for (auto it = gene_collection.genes.begin();
//      it != gene_collection.genes.end(); ++it) {
////        make_dir(folder + ToString(it->first));
//    const GenePositions& gene_poss = it->second.gene_positions;
//
//    //todo improve later
//    Sequence total_gene_sequence;
//    for (GenomeId genome_id : key_set(gene_collection.genomes)) {
//      const Sequence& genome = get(gene_collection.genomes, genome_id).sequence;
//      for (Pos pos : get_all(gene_poss, genome_id)) {
//        total_gene_sequence = total_gene_sequence + genome.Subseq(pos.first.start_pos, pos.first.end_pos);
//      }
//    }
//    WriteComponentsAlongSequence(gp, folder + ToString(it->first) + "/",
//                                 100000, 50, total_gene_sequence, coloring);
//  }
//}

}
