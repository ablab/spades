#pragma once

#include "transitions.hpp"
#include "contracted_graph.hpp"
#include "cluster_storage_builder.hpp"
#include "cluster_storage_analyzer.hpp"

namespace contracted_graph {
    struct ContractedVertexData {
        size_t indegree_;
        size_t outdegree_;
        size_t capacity_;
    };

    class VertexDataStats: public read_cloud_statistics::Statistic {
        std::unordered_map<VertexId, ContractedVertexData> data_;

     public:
        VertexDataStats(): read_cloud_statistics::Statistic("contracted_vertex_stats"), data_() {}
        VertexDataStats(const VertexDataStats& other) = default;
        void Insert(const VertexId& vertex, size_t indegree, size_t outdegree, size_t capacity) {
            data_[vertex].indegree_ = indegree;
            data_[vertex].outdegree_ = outdegree;
            data_[vertex].capacity_ = capacity;
        }

        void Serialize(ofstream& fout) override {
            for (const auto& entry: data_) {
                fout << entry.first.int_id() << " " << entry.second.indegree_ << " "
                     << entry.second.outdegree_ << " " << entry.second.capacity_ << std::endl;
            }
        }
    };

    struct VertexClusterData {
      std::string name_;
      std::size_t clusters_;
      std::size_t correct_clusters_;
      std::map<transitions::Transition, size_t> transition_to_clusters_;
      std::map<transitions::Transition, size_t> correct_transition_to_clusters_;
      bool has_correct_;

      VertexClusterData(const string &name = "", bool has_correct = true)
          : name_(name), clusters_(0), correct_clusters_(0), transition_to_clusters_(),
            correct_transition_to_clusters_(), has_correct_(has_correct) {}

      void AddTransitionWithCluster(const transitions::Transition& transition, size_t clusters) {
          transition_to_clusters_[transition] = clusters;
      }

      void AddCorrectTransitionWithCluster(const transitions::Transition& transition, size_t clusters) {
          correct_transition_to_clusters_[transition] = clusters;
      }

      void Serialize(ofstream& fout, const string& sep) const {
          fout << clusters_;
          if (clusters_ != 0) {
              fout << "=";
              PrintMap(fout, transition_to_clusters_);
          }
          if (has_correct_) {
              fout << sep;
              fout << correct_clusters_;
              if (correct_clusters_ != 0) {
                  fout << "=";
                  PrintMap(fout, correct_transition_to_clusters_);
              }
          }
          fout << sep;
      }

      void PrintMap(ofstream& fout, const std::map<transitions::Transition, size_t>& transition_map) const {
          vector<size_t> weights;
          string buffer;
          for (const auto &transition_entry: transition_map) {
              weights.push_back(transition_entry.second);
          }
          for (const auto& weight: weights) {
              buffer += std::to_string(weight) + "+";
          }
          if (buffer.size() != 0) {
              buffer.pop_back();
          }
          fout << buffer;
      }
    };

    struct TransitionVertexData {
      vector<std::pair<size_t, size_t>> transitions_;
      std::set<EdgeId> incoming_set_;
      std::set<EdgeId> outcoming_set_;
      std::map<EdgeId, size_t> incoming_edge_to_pos_;
      std::map<EdgeId, size_t> outcoming_edge_to_pos_;
      string name_;

      TransitionVertexData(const set<EdgeId> &incoming_set_, const set<EdgeId> &outcoming_set_, const string& name) :
          transitions_(), incoming_set_(incoming_set_), outcoming_set_(outcoming_set_),
          incoming_edge_to_pos_(), outcoming_edge_to_pos_(), name_(name) {
          size_t counter = 1;
          for (const auto& edge: incoming_set_) {
              incoming_edge_to_pos_[edge] = counter;
              ++counter;
          }
          counter = 1;
          for (const auto& edge: outcoming_set_) {
              outcoming_edge_to_pos_[edge] = counter;
              ++counter;
          }
      }

      TransitionVertexData(const TransitionVertexData& other) = default;

      void InsertTransition(const EdgeId& incoming, const EdgeId& outcoming) {
          size_t incoming_index = incoming_edge_to_pos_.at(incoming);
          size_t outcoming_index = outcoming_edge_to_pos_.at(outcoming);
          transitions_.push_back({incoming_index, outcoming_index});
      }

      void Print(ofstream& fout, const string& separator) const {
          fout << transitions_.size();
          if (transitions_.size() != 0) {
              fout << "(";
              string buffer;
              for (const auto &transition: transitions_) {
                  buffer += "I";
                  buffer += std::to_string(transition.first) + "-" + "O" + std::to_string(transition.second) + ",";
              }
              buffer.pop_back();
              fout << buffer;
              fout << ")";
          }
          fout << separator;
      }

      string GetName() const {
          return name_;
      }

    };

    struct DetailedVertexData {
        vector<size_t> lengths_;
        vector<double> coverages_;
        vector<size_t> barcode_coverages_;
        size_t incoming_;
        size_t outcoming_;
        size_t capacity_;
        vector<TransitionVertexData> vertex_to_transition_data_;
        vector<VertexClusterData> vertex_to_cluster_data_;
        VertexClusterData vertex_to_nonpath_data_;
    };

    class TransitionStatistics: public read_cloud_statistics::Statistic {
        std::map<VertexId, DetailedVertexData> data_;

     public:
        typedef std::map<VertexId, DetailedVertexData>::const_iterator const_iterator;
        TransitionStatistics()
            : Statistic("transition_statistics"), data_() {}

        void Serialize(ofstream& fout) override {
            string sep = "\t";
            fout << "Incoming" << sep << "Outcoming" << sep << "Capacity" << sep;
            fout << GetTransitionColumnNames(sep);
            fout << "Length" << sep << "Coverages" << sep << "Barcode coverages" << sep;
            fout << GetClusterColumnNames(sep);
            fout << GetNonpathColumnName(sep);
            fout << std::endl;

            for (const auto& entry: data_) {
                if (entry.second.incoming_ == 2 and entry.second.outcoming_ == 2) {
                    fout << entry.second.incoming_ << sep << entry.second.outcoming_ << sep << entry.second.capacity_
                         << sep;
                    PrintTransitionData(entry.second.vertex_to_transition_data_, fout, sep);
                    PrintLengths(entry.second.lengths_, fout, sep);
                    PrintCoverages(entry.second.coverages_, fout, sep);
                    PrintVectorWithSeparator(entry.second.barcode_coverages_, fout, sep);
                    for (const auto &cluster_data: entry.second.vertex_to_cluster_data_) {
                        cluster_data.Serialize(fout, sep);
                    }
                    entry.second.vertex_to_nonpath_data_.Serialize(fout, sep);
                    fout << std::endl;
                }
            }
        }

        string GetTransitionColumnNames(const string& sep) {
            VERIFY(data_.size() != 0);
            DetailedVertexData example = (*(data_.begin())).second;
            string result;
            for (const auto& transition_data: example.vertex_to_transition_data_) {
                result += transition_data.GetName() + " transitions" + sep;
            }
            return result;
        }

        string GetClusterColumnNames(const string& sep) {
            VERIFY(data_.size() != 0);
            DetailedVertexData example = (*(data_.begin())).second;
            string result;
            for (const auto& cluster_data: example.vertex_to_cluster_data_) {
                result += cluster_data.name_;
                result += sep;
                result += "Correct " + cluster_data.name_;
                result += sep;
            }
            return result;
        }

        string GetNonpathColumnName(const string& sep) {
            VERIFY(data_.size() != 0);
            DetailedVertexData example = (*(data_.begin())).second;
            string result;
            VertexClusterData cluster_data = example.vertex_to_nonpath_data_;
            result += cluster_data.name_;
            result += sep;
            return result;
        }

        void PrintTransitionData(const vector<TransitionVertexData>& transition_data, ofstream& fout, const string separator) {
            for (const auto& transition_storage: transition_data) {
                transition_storage.Print(fout, separator);
            }
        }

        void PrintLengths(const vector<size_t>& lengths, ofstream& fout, const string separator) {
            vector<size_t> new_lengths;
            for (const auto& length: lengths) {
                new_lengths.push_back(length / 1000);
            }
            PrintVectorWithSeparator(new_lengths, fout, separator);
        }

        void PrintCoverages(const vector<double>& coverages, ofstream& fout, const string separator) {
            vector<size_t> new_coverages;
            for (const auto& coverage: coverages) {
                new_coverages.push_back(static_cast<size_t>(coverage));
            }
            PrintVectorWithSeparator(new_coverages, fout, separator);
        }

        template<class T>
        void PrintVectorWithSeparator(const vector<T>& vec, ofstream& fout, const string external_separator,
                                      const string internal_separator = ",") {
            std::string buffer;
            for (const auto& el: vec) {
                buffer += std::to_string(el) + internal_separator;
            }
            buffer.pop_back();
            buffer += external_separator;
            fout << buffer;
        }

        void Insert(const VertexId& vertex, const DetailedVertexData& statistics) {
            data_[vertex] = statistics;
        }

        const_iterator begin() const {
            return data_.begin();
        }

        const_iterator end() const {
            return data_.end();
        }
    };

    struct TransitionResolutionStats : public read_cloud_statistics::Statistic {
      size_t reference_resolved_transitions_;
      size_t contig_resolved_transitions_;
      size_t cloud_resolved_transitions_;
      size_t resolved_contig_vertices_;
      size_t resolved_cloud_vertices_;
      size_t non_reference_contig_transitions_;
      size_t non_reference_cloud_transitions_;

      TransitionResolutionStats()
          : Statistic("transition_resolution_statistics"),
            reference_resolved_transitions_(0),
            contig_resolved_transitions_(0),
            cloud_resolved_transitions_(0),
            resolved_contig_vertices_(0),
            resolved_cloud_vertices_(0),
            non_reference_contig_transitions_(0),
            non_reference_cloud_transitions_(0) {}

      void Serialize(ofstream& fout) override {
          fout << "Reference resolved transitions: " << reference_resolved_transitions_ << "\n";
          fout << "Contig resolved transitions: " << contig_resolved_transitions_ << "\n";
          fout << "Cloud resolved transitions: " << cloud_resolved_transitions_ << "\n";
          fout << "Resolved contig vertices" << resolved_contig_vertices_ << "\n";
          fout << "Resolved cloud vertices: " << resolved_cloud_vertices_ << "\n";
          fout << "False contig transitions: " << non_reference_contig_transitions_ << "\n";
          fout << "False cloud transitions: " << non_reference_cloud_transitions_ << std::endl;
      }
    };

    class ContractedGraphAnalyzer: public read_cloud_statistics::StatisticProcessor {
        const Graph& graph_;
        const FrameBarcodeIndexInfoExtractor& barcode_extractor_;
        const cluster_statistics::PathClusterStorage& path_cluster_storage_;
        contracted_graph::ContractedGraph contracted_graph_;
        const unordered_map<string, transitions::ContigTransitionStorage>& name_to_transition_storage_;
        const transitions::ContigTransitionStorage& reference_transition_storage_;
        const cluster_statistics::ClusterStorage& cluster_storage_;
        const size_t min_read_threshold_;
     public:
        ContractedGraphAnalyzer(const Graph& graph, const FrameBarcodeIndexInfoExtractor& barcode_extractor,
                                const cluster_statistics::PathClusterStorage &path_cluster_storage_,
                                const ContractedGraph &contracted_graph_,
                                const unordered_map<string, transitions::ContigTransitionStorage>& name_to_storage,
                                const transitions::ContigTransitionStorage& reference_transition_storage,
                                const cluster_statistics::ClusterStorage& cluster_storage,
                                size_t min_read_threshold) :
            StatisticProcessor("transition_statistics_extractor"), graph_(graph), barcode_extractor_(barcode_extractor),
            path_cluster_storage_(path_cluster_storage_), contracted_graph_(contracted_graph_),
            name_to_transition_storage_(name_to_storage),
            reference_transition_storage_(reference_transition_storage),
            cluster_storage_(cluster_storage),
            min_read_threshold_(min_read_threshold) {}

        void FillStatistics() override {
            auto transition_statistics_ptr = make_shared<TransitionStatistics>(GetTransitionStatistics(min_read_threshold_));
            auto vertex_statistics_ptr = make_shared<VertexDataStats>(GetVertexDataStats());
            AddStatistic(transition_statistics_ptr);
            AddStatistic(vertex_statistics_ptr);
        }

        DECL_LOGGER("ContractedGraphAnalyzer");

     private:
        VertexDataStats GetVertexDataStats() {
            VertexDataStats stats;
            for (const auto& entry: contracted_graph_) {
                VertexId vertex = entry.first;
                size_t outcoming = contracted_graph_.GetOutcoming(vertex).size();
                size_t incoming = contracted_graph_.GetIncoming(vertex).size();
                size_t capacity = contracted_graph_.GetCapacity(vertex);
                stats.Insert(vertex, incoming, outcoming, capacity);
            }
            return stats;
        }

        TransitionStatistics GetTransitionStatistics(size_t min_read_threshold) {
            TransitionStatistics stats;
            auto cluster_storage_map = BuildTransitionClusterStorages(path_cluster_storage_);
            auto edge_cluster_storage = BuildEdgeClusterStorage(cluster_storage_);
            INFO("Getting transition statistics...")
            for (const auto& entry: contracted_graph_) {
                VertexId vertex = entry.first;
                DEBUG("Vertex: " << vertex.int_id());
                DetailedVertexData vertex_data;
                vector<EdgeId> incoming = contracted_graph_.GetIncoming(vertex);
                vector<EdgeId> outcoming = contracted_graph_.GetOutcoming(vertex);
                vertex_data.incoming_ = incoming.size();
                vertex_data.outcoming_ = outcoming.size();
                vertex_data.capacity_ = contracted_graph_.GetCapacity(vertex);

                vector<VertexClusterData> vertex_to_cluster_storages = {VertexClusterData("Path clusters assigned to the vertex (2 edges)"),
                                                                        VertexClusterData("Path clusters assigned to the vertex (3 edges)"),
                                                                        VertexClusterData("Path clusters assigned to the vertex (>3 edges)")};

                VertexClusterData nonpath_cluster_storage("Non-path clusters assigned to the vertex", false);

                std::set<EdgeId> incoming_set(incoming.begin(), incoming.end());
                std::set<EdgeId> outcoming_set(outcoming.begin(), outcoming.end());

                vector<TransitionVertexData> vertex_to_transition_storages;
                for (const auto& storage_entry: name_to_transition_storage_) {
                    auto transition_storage = TransitionVertexData(incoming_set, outcoming_set, storage_entry.first);
                    vertex_to_transition_storages.push_back(transition_storage);
                }

                //fixme why?
                const size_t coverage_read_threshold = 3;
                for (const auto& edge: incoming_set) {
                    vertex_data.lengths_.push_back(graph_.length(edge));
                    vertex_data.coverages_.push_back(graph_.coverage(edge));
                    vertex_data.barcode_coverages_.push_back(GetBarcodeCoverage(edge, coverage_read_threshold));
                }
                for (const auto& edge: outcoming_set) {
                    vertex_data.lengths_.push_back(graph_.length(edge));
                    vertex_data.coverages_.push_back(graph_.coverage(edge));
                    vertex_data.barcode_coverages_.push_back(GetBarcodeCoverage(edge, coverage_read_threshold));
                }
                for (const auto& in_edge: incoming_set) {
                    for (const auto& out_edge: outcoming_set) {
                        DEBUG("Incoming edge: " << in_edge.int_id());
                        DEBUG("Outcoming edge: " << out_edge.int_id());
                        vector<EdgeId> edges = {in_edge, out_edge};
                        cluster_statistics::SimplePath path(edges);
                        DEBUG("Number of clusters");
                        transitions::Transition transition(in_edge, out_edge);
                        UpdateTransitionStorages(vertex_to_transition_storages, in_edge, out_edge);

                        for (auto& vertex_cluster_data: vertex_to_cluster_storages) {
                            UpdateTransitionToClusters(vertex_cluster_data, transition,
                                                       cluster_storage_map[vertex_cluster_data.name_], min_read_threshold);
                        }
                        UpdateNonPathClusters(nonpath_cluster_storage, transition, edge_cluster_storage);
                    }
                }
                vertex_data.vertex_to_cluster_data_ = std::move(vertex_to_cluster_storages);
                vertex_data.vertex_to_nonpath_data_ = std::move(nonpath_cluster_storage);
                vertex_data.vertex_to_transition_data_ = vertex_to_transition_storages;
                stats.Insert(vertex, vertex_data);
            }
            return stats;
        }

        size_t GetBarcodeCoverage(const EdgeId& edge, const size_t read_threshold) {
            vector<BarcodeId> barcodes = barcode_extractor_.GetBarcodes(edge);
            size_t barcode_coverage = (size_t) std::count_if(barcodes.begin(), barcodes.end(),
                                                             [read_threshold, &edge, this](const BarcodeId& barcode) {
                                                               return barcode_extractor_.GetNumberOfReads(edge, barcode) > read_threshold;
                                                             });
            return barcode_coverage;
        }

        void UpdateTransitionStorages(vector<TransitionVertexData>& transition_storages,
                                      const EdgeId& in_edge, const EdgeId& out_edge) {
            DEBUG("Transition to storages")
            for (auto& storage: transition_storages) {
                string name = storage.GetName();
                if (name_to_transition_storage_.at(name).CheckTransition(in_edge, out_edge)) {
                    storage.InsertTransition(in_edge, out_edge);
                }
            }
        }

        void UpdateTransitionToClusters(VertexClusterData& vertex_cluster_data, const transitions::Transition& transition,
                                        const cluster_statistics::TransitionClusterStorage& transition_cluster_storage,
                                        size_t min_read_threshold) {
            size_t clusters = GetNumberOfClusters(transition, transition_cluster_storage, min_read_threshold);
            vertex_cluster_data.clusters_ += clusters;
            vertex_cluster_data.transition_to_clusters_[transition] = clusters;
            DEBUG("Checking transition");
            if (reference_transition_storage_.CheckTransition(transition)) {
                vertex_cluster_data.correct_clusters_ += clusters;
                vertex_cluster_data.correct_transition_to_clusters_[transition] = clusters;
            }
        }

        void UpdateNonPathClusters(VertexClusterData& vertex_cluster_data, const transitions::Transition& transition,
                                   const cluster_statistics::EdgeClusterStorage& edge_cluster_storage) {
            DEBUG("Non path clusters")
            const EdgeId first = transition.first_;
            const EdgeId second = transition.second_;
            if (edge_cluster_storage.HasKey(first) and edge_cluster_storage.HasKey(second)) {
                vector<size_t> first_cluster_ids;
                vector<size_t> second_cluster_ids;
                std::for_each(edge_cluster_storage.begin(first), edge_cluster_storage.end(first),
                              [&first_cluster_ids](const cluster_statistics::Cluster &cluster) {
                                first_cluster_ids.push_back(cluster.GetId());
                              });
                std::for_each(edge_cluster_storage.begin(second), edge_cluster_storage.end(second),
                              [&second_cluster_ids](const cluster_statistics::Cluster &cluster) {
                                second_cluster_ids.push_back(cluster.GetId());
                              });
                std::sort(first_cluster_ids.begin(), first_cluster_ids.end());
                std::sort(second_cluster_ids.begin(), second_cluster_ids.end());
                vector<size_t> intersection;
                std::set_intersection(first_cluster_ids.begin(), first_cluster_ids.end(),
                                      second_cluster_ids.begin(), second_cluster_ids.end(),
                                      std::back_inserter(intersection));
                size_t clusters = intersection.size();
                vertex_cluster_data.clusters_ += clusters;
                vertex_cluster_data.transition_to_clusters_[transition] = clusters;
            }
            DEBUG("Finished building non path clusters")
        }

        std::unordered_map<string, cluster_statistics::TransitionClusterStorage>
        BuildTransitionClusterStorages(const cluster_statistics::PathClusterStorage path_cluster_storage) {
            std::unordered_map<string, cluster_statistics::TransitionClusterStorage> result;
            cluster_statistics::PredicateTransitionClusterStorageBuilder builder;
            cluster_statistics::TwoEdgePathPredicate two_edge_path_predicate;
            cluster_statistics::ThreeEdgePathPredicate three_edge_path_predicate;
            cluster_statistics::ManyEdgePathPredicate many_edge_path_predicate;
            result.insert({"Path clusters assigned to the vertex (2 edges)",
                           builder.BuildTransitionClusterStorage(path_cluster_storage, two_edge_path_predicate)});
            result.insert({"Path clusters assigned to the vertex (3 edges)",
                           builder.BuildTransitionClusterStorage(path_cluster_storage, three_edge_path_predicate)});
            result.insert({"Path clusters assigned to the vertex (>3 edges)",
                           builder.BuildTransitionClusterStorage(path_cluster_storage, many_edge_path_predicate)});
            return result;
        };

        cluster_statistics::EdgeClusterStorage BuildEdgeClusterStorage(const cluster_statistics::ClusterStorage& cluster_storage) {
            INFO("Edge cluster storage")
            cluster_statistics::EdgeClusterStorageBuilder edge_cluster_storage_builder;
            return edge_cluster_storage_builder.BuildEdgeClusterStorage(cluster_storage, min_read_threshold_);
        }

        size_t GetNumberOfClusters(const transitions::Transition& transition,
                                   const cluster_statistics::TransitionClusterStorage& transition_cluster_storage,
                                   size_t min_read_threshold) const {
            if (transition_cluster_storage.HasKey(transition)) {
                size_t clusters = (size_t) std::count_if(transition_cluster_storage.begin(transition),
                                                         transition_cluster_storage.end(transition),
                                                         [min_read_threshold](const cluster_statistics::Cluster &cluster) {
                                                           return cluster.GetReads() > min_read_threshold;
                                                         });
                return clusters;
            }
            return 0;
        }

//        TransitionResolutionStats GetTransitionResolutionStats(const TransitionStatistics& transition_statistics) {
//            TransitionResolutionStats result;
//            for (const auto& entry: transition_statistics) {
//                const vector<TransitionVertexData>& transition_storages = entry.second.vertex_to_transition_data_;
//                auto reference_transitions = FindTransitionDataByName(transition_storages, "Reference").transitions_;
//                auto contig_transitions = FindTransitionDataByName(transition_storages, "Contig").transitions_;
//                auto cloud_contig_transitions = FindTransitionDataByName(transition_storages, "Read cloud contig").transitions_;
//
//            }
//        }

//        TransitionVertexData FindTransitionDataByName(const vector<TransitionVertexData>& transition_storages,
//                                                      const string& name) {
//            auto it = std::find_if(transition_storages.begin(), transition_storages.end(), [&name](TransitionVertexData& transition_data) {
//              return transition_data.GetName() == name;
//            });
//            VERIFY(it != transition_storages.end());
//            return *it;
//        }
    };
}