//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "alignment_mixer.hpp"
#include "common/assembly_graph/paths/bidirectional_path_io/io_support.hpp"
#include "common/io/binary/graph_pack.hpp"
#include "helpers/aligner_output_reader.hpp"
#include "helpers/common.hpp"
#include "helpers/replacer.hpp"
#include "helpers/replace_info_writer.hpp"

#include "utils/filesystem/file_opener.hpp"
#include "utils/filesystem/path_helper.hpp"

#include <string>
#include <algorithm>
#include <unordered_map>

namespace helpers {

enum class GAFColumns : size_t {
    contig_name            = 0,
    contig_len             = 1,
    contig_start_pos       = 2, // (0-based; closed)
    contig_end_pos         = 3, // (0-based; open)
    strand                 = 4,
    path                   = 5,
    path_len               = 6,
    path_start_pos         = 7,
    path_end_pos           = 8,
    matches                = 9,
    alignment_block_length = 10,
    mapping_quality        = 11,
    nm_tag                 = 12,
    as_tag                 = 13,
    dv_tag                 = 14,
    id_tag                 = 15,
    sg_tag                 = 16,
    TOTAL_COLUMNS_SIZE     = 17
};

template<GAFColumns el>
struct type_getter<GAFColumns, el> {
    using type = 
        std::conditional_t<
            el == GAFColumns::strand,
            char,
        std::conditional_t<
            el == GAFColumns::contig_len             ||
            el == GAFColumns::contig_start_pos       ||
            el == GAFColumns::contig_end_pos         ||
            el == GAFColumns::path_len               ||
            el == GAFColumns::path_start_pos         ||
            el == GAFColumns::path_end_pos           ||
            el == GAFColumns::matches                ||
            el == GAFColumns::alignment_block_length ||
            el == GAFColumns::mapping_quality,
            unsigned long long,
            std::string>>;
};

template<GAFColumns ... columns>
Records<GAFColumns, columns ...> ReadGraphAlignerOutput(std::string const & file) {
    auto inp = fs::open_file(file, std::ios_base::in, std::ios_base::badbit);
    auto f = FilterType<GAFColumns, columns ...>([](auto const & record){ return record.template Get<GAFColumns::strand>() == '+'; });
    Records<GAFColumns, columns ...> records;
    RecordPusher<GAFColumns, columns ...> pusher(records, f);
    std::string line;
    size_t total_lines = 0;
    size_t accepted_lines = 0;
    while (GetNextNonemptyLine(inp, line)) {
        accepted_lines += pusher.Push(line);
        ++total_lines;
    }
    INFO("Total line read: " << total_lines);
    INFO("Accepted lines: " << accepted_lines);
    return records;
}

} // namespace helpers

namespace alignment_mixer {

using namespace debruijn_graph;
using namespace path_extend;
using namespace std;
using namespace helpers;

struct gcfg {
    std::string contigs_file;
    std::string graph_file;
    std::string alignments;
    std::string output_dir;
    size_t k;
} cfg;

clipp::group GetCLI() {
  using namespace clipp;

  auto cli = (
      cfg.contigs_file << value("contigs file"),
      cfg.graph_file << value("graph.gfa"),
      cfg.alignments << value("alignments.gaf"),
      cfg.output_dir << value("output dir"),
      (required("-k") & integer("int", cfg.k)) % "k-mer length to use"
  );

  return cli;
}

namespace {

EdgeId getEdgeId(std::string const & s, size_t & pos) {
    size_t len = 0;
    while (pos + len < s.size() && isdigit(s[pos + len]))
        ++len;
    auto ans = std::stoull(s.substr(pos, len));
    pos += len;
    return ans;
}

std::string MakeSequence(debruijn_graph::Graph const & graph, std::string const & path, size_t from, size_t to, size_t len) {
    size_t pos = 0;
    std::vector<EdgeId> edges;
    while (pos + 1 < path.size()) {
        char direction = path[pos];
        auto edge_id = getEdgeId(path, ++pos);
        switch (direction) {
        case '>': edges.push_back(edge_id); break;
        case '<': edges.push_back(graph.conjugate(edge_id)); break;
        default: throw "wrong path format";
        }
    }
    ScaffoldSequenceMaker seq_maker(graph);
    auto query_seq = seq_maker.MakeSequence(*BidirectionalPath::create(graph, edges));
    VERIFY(len == query_seq.size());
    return query_seq.substr(from, to - from);
}


template <GAFColumns ... columns>
std::string MakeSequence(debruijn_graph::Graph const & graph, Records<GAFColumns, columns ...> const & records, size_t id) {
    auto const & start_pos = records[id].template Get<GAFColumns::path_start_pos>();
    auto const & end_pos = records[id].template Get<GAFColumns::path_end_pos>();
    auto const & len = records[id].template Get<GAFColumns::path_len>();
    auto const & path = records[id].template Get<GAFColumns::path>();
    return MakeSequence(graph, path, start_pos, end_pos, len);
}

std::string TagContent(std::string const & tag) {
    return tag.substr(5);
}

/*
bool Contains(std::vector<size_t> const & vs, size_t v) {
    return std::find(vs.begin(), vs.end(), v) != vs.end();
}

bool Contains(std::vector<size_t> const & vs, std::vector<size_t> const & elements) {
    for (auto const & e : elements) {
        if (!Contains(vs, e))
            return false;
    }
    return true;
}

template <GAFColumns ... columns>
void Print(std::vector<size_t> const & vs, Records<GAFColumns, columns ...> const & records) {
    for (auto const & id : vs) {
        cout << id << ": "
                << records[id].template Get<GAFColumns::contig_start_pos>() << " -> "
                << records[id].template Get<GAFColumns::contig_end_pos>() << '\n';
    }
}

template <GAFColumns ... columns>
void MakeClusters(Records<GAFColumns, columns ...> const & records) {
    std::unordered_map<std::string, std::unordered_map<size_t, std::vector<size_t>>> starts_from_same_pos;
    std::unordered_map<std::string, std::unordered_map<size_t, std::vector<size_t>>> ends_on_same_pos;

    auto Check = [&](std::string const & name, auto const & records, std::vector<size_t> const & cluster, std::vector<size_t> const & first_set, std::vector<size_t> const & second_set) {
        auto Printer = [&, is_bad_sets = false](auto const & id, auto const & another_set) mutable {
            if (Contains(cluster, another_set))
                return;
            if (!is_bad_sets) {
                is_bad_sets = true;
                cout << name << ":\n";
                Print(first_set, records);
                cout << "---------------\n";
                Print(second_set, records);
                cout << "===============\n";
            }
            cout << "id: " << id << '\n';
            Print(another_set, records);
            cout << "---------------\n";
        };

        for (auto const & id : first_set) {
            if (Contains(second_set, id))
                continue;
            auto end_pos = records[id].template Get<GAFColumns::contig_end_pos>();
            auto const & another_set = ends_on_same_pos[name][end_pos];
            Printer(id, another_set);
        }
        for (auto const & id : second_set) {
            if (Contains(first_set, id))
                continue;
            auto start_pos = records[id].template Get<GAFColumns::contig_start_pos>();
            auto const & another_set = starts_from_same_pos[name][start_pos];
            Printer(id, another_set);
        }
    };

    for (size_t i = 0; i < records.size(); ++i) {
        auto const & record = records[i];
        auto const & contig_name = record.template Get<GAFColumns::contig_name>();
        auto const & contig_start_pos = record.template Get<GAFColumns::contig_start_pos>();
        auto const & contig_end_pos = record.template Get<GAFColumns::contig_end_pos>();
        starts_from_same_pos[contig_name][contig_start_pos].push_back(i);
        ends_on_same_pos[contig_name][contig_end_pos].push_back(i);
    }

    std::unordered_map<std::string, std::vector<std::vector<size_t>>> clusters;
    for (auto const & starts : starts_from_same_pos) {
        for (auto const & pos_and_indices : starts.second) {
            auto const & first_set = pos_and_indices.second;
            auto end_pos = records[first_set.front()].template Get<GAFColumns::contig_end_pos>();
            auto const & second_set = ends_on_same_pos[starts.first][end_pos];
            if (first_set == second_set) {
                clusters[starts.first].push_back(first_set);
                continue;
            }
            std::vector<size_t> cluster;
            std::set_union(first_set.begin(), first_set.end(), second_set.begin(), second_set.end(), std::back_inserter(cluster));
            Check(starts.first, records, cluster, first_set, second_set);
            clusters[starts.first].push_back(std::move(cluster));
        }
    }
}
*/

using NamedGroups = std::unordered_map<std::string, std::vector<size_t>>;

struct ScanlineElement {
    size_t id;
    size_t pos;
    bool ending;

    bool operator < (ScanlineElement const & other) const noexcept {
        if (pos != other.pos)
            return pos < other.pos;
        if (ending != other.ending)
            return ending > other.ending;
        return id < other.id;
    }
};

using Cluster = std::vector<size_t>;
using Clusters = std::vector<Cluster>;


template <GAFColumns ... columns>
void Print(std::vector<ScanlineElement> const & vs, Records<GAFColumns, columns ...> const & records) {
    for (auto const & el : vs) {
        cout << el.id << "( " << el.ending << " ): "
                << records[el.id].template Get<GAFColumns::contig_start_pos>() << " -> "
                << records[el.id].template Get<GAFColumns::contig_end_pos>() << '\n';
    }
    cout << "--------------\n";
}

template <GAFColumns ... columns>
Clusters MakeClusters(Records<GAFColumns, columns ...> const & records, std::vector<size_t> const & indices) {
    auto StartPos = [&records](auto const & id) { return records[id].template Get<GAFColumns::contig_start_pos>(); };
    auto EndPos = [&records](auto const & id) { return records[id].template Get<GAFColumns::contig_end_pos>(); };
    std::vector<ScanlineElement> elements;
    elements.reserve(indices.size() * 2);

    for (auto const & id : indices) {
        elements.push_back({id, StartPos(id), false});
        elements.push_back({id, EndPos(id), true});
    }

    std::sort(elements.begin(), elements.end());

    Clusters clusters;
    size_t nesting = 0;
    Cluster cluster;
    bool bad_clusters = false;
    bool may_increase = true;
    for (auto const & el : elements) {
        if (el.ending) {
            may_increase = false;
            --nesting;
            if (nesting == 0) {
                may_increase = true;
                clusters.push_back(std::move(cluster));
                cluster.clear();
            }
        } else {
            if (!may_increase) {
                bad_clusters = true;
                // cout << records[el.id].template Get<GAFColumns::contig_name>() << '\n';
            }
            ++nesting;
            cluster.push_back(el.id);
        }
    }

    // if (bad_clusters)
    //     Print(elements, records);
    return clusters;
}

template <GAFColumns ... columns>
NamedGroups GroupByName(Records<GAFColumns, columns ...> const & records) {
    NamedGroups groups;
    for (size_t i = 0; i < records.size(); ++i) {
        auto const & name = records[i].template Get<GAFColumns::contig_name>();
        groups[name].push_back(i);
    }
    return groups;
}

template <GAFColumns ... columns>
ReplaceInfo Mix(debruijn_graph::Graph const & graph, Records<GAFColumns, columns ...> const & records, Cluster const & cluster) {
    auto StartPos = [&records](auto const & id) { return records[id].template Get<GAFColumns::contig_start_pos>(); };
    auto EndPos = [&records](auto const & id) { return records[id].template Get<GAFColumns::contig_end_pos>(); };
    auto GetScore = [&records](size_t const & id) { return stod(TagContent(records[id].template Get<GAFColumns::id_tag>())); };
    size_t best_path_id = cluster.front();
    double best_score = GetScore(best_path_id);
    for (size_t i = 1; i < cluster.size(); ++i) {
        auto id = cluster[i];
        auto score = GetScore(id);
        if (score > best_score) {
            best_score = score;
            best_path_id = id;
        }
    }

    return ReplaceInfo(MakeSequence(graph, records, best_path_id), StartPos(best_path_id), EndPos(best_path_id));
    // auto best_path_seq = ;
    
}

void WriteWithWidth(std::ostream & out, std::string const & seq, size_t width = 50) {
    for (size_t pos = 0; pos < seq.size(); pos += width)
        out << seq.substr(pos, width) << '\n';
}

#define WISHED_COLIMNS GAFColumns::contig_name, GAFColumns::contig_start_pos, GAFColumns::contig_end_pos,\
    GAFColumns::path, GAFColumns::path_len, GAFColumns::path_start_pos, GAFColumns::path_end_pos, GAFColumns::id_tag, GAFColumns::sg_tag, GAFColumns::strand

constexpr char BASE_NAME[] = "graph_pack";

} //namespace

int main() {
    START_BANNER("SPAdes standalone aligner mixer");
    ReplaceInfoWriter::SetStream(fs::append_path(cfg.output_dir, "replace_info_dump.dump"));
    auto k = cfg.k;
    auto contigs = ReadContigs(cfg.contigs_file);
    auto paths = helpers::ReadGraphAlignerOutput<WISHED_COLIMNS>(cfg.alignments);
    auto groups = GroupByName(paths);
    std::unordered_map<std::string, Clusters> namedClusters;
    for (auto const & group : groups)
        namedClusters[group.first] = MakeClusters(paths, group.second);

    CHECK_FATAL_ERROR(runtime_k::MIN_K <= k, "k-mer size " << k << " is too low");
    CHECK_FATAL_ERROR(k < runtime_k::MAX_K, "k-mer size " << k << " is too high, recompile with larger SPADES_MAX_K option");
    CHECK_FATAL_ERROR(k & 1, "k-mer size must be odd");

    fs::make_dir(cfg.output_dir);
    debruijn_graph::GraphPack gp(k, cfg.output_dir, 0);
    auto p = fs::append_path(cfg.graph_file, BASE_NAME);
    io::binary::FullPackIO().Load(p, gp);
    auto const & graph = gp.get<Graph>();

    for (auto & contig : contigs) {
        auto clusters = namedClusters.find(contig.name);
        if (clusters == namedClusters.end())
            continue;
        std::list<ReplaceInfo> infos;
        for (auto const & cluster : clusters->second)
            infos.push_back(Mix(graph, paths, cluster));
        contig.seq = ReplaceAndDump(contig.seq, std::move(infos), contig.name);
        contig.corrected = true;
    }

    ofstream contigs_output(fs::append_path(cfg.output_dir, "corrected.fasta"));
    VERIFY(contigs_output.is_open());
    for (auto const & contig : contigs) {
        if (contig.seq.empty())
            continue;
        contigs_output << '>' << contig.name << " len=" << contig.seq.size() << " is_corrected=" << contig.corrected << '\n';
        WriteWithWidth(contigs_output, contig.seq);
    }

    INFO("SPAdes standalone aligner mixer finished");
    return 0;
}

} // namespace alignment_mixer
