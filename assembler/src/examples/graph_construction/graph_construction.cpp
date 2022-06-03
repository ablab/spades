#include "graph_construction.hpp"

#include "common/utils/filesystem/temporary.hpp"
#include "common/io/reads/coverage_filtering_read_wrapper.hpp" // io::CovFilteringWrap
#include "common/io/reads/read_stream_vector.hpp" // io::ReadStreamList
#include "common/io/reads/single_read.hpp" // io::SingleReadSeq
#include "common/assembly_graph/construction/debruijn_graph_constructor.hpp"
#include "common/kmer_index/extension_index/kmer_extension_index.hpp"
#include "common/kmer_index/extension_index/kmer_extension_index_builder.hpp"
#include "common/kmer_index/kmer_mph/kmer_index_builder.hpp"
#include "common/configs/config_struct.hpp"
#include "common/kmer_index/kmer_mph/kmer_splitters.hpp"
#include "examples/graph_io/gfa_io.hpp"

#include "common/utils/logger/log_writers.hpp"

#include <random>
#include <sstream>

namespace spades_example {

    class StreamWrapper {
    public:
        StreamWrapper(const std::vector<std::string>& reads)
            : reads_(reads), idx_(0), is_open_(1) {}
        bool is_open() { return is_open_; }
        bool eof() { return idx_ == reads_.size(); }
        void close() { is_open_ = 0; }
        void reset() { idx_ = 0; }
        StreamWrapper& operator>>(io::SingleReadSeq& read) {
            read = io::SingleReadSeq(Sequence(reads_[idx_++]));
            return *this;
        }
    private:
        std::vector<std::string> reads_;
        size_t idx_;
        bool is_open_;
    };

void CreateRandomGraph(debruijn_graph::Graph& g, size_t genome_length) {
    INFO("start creation");
    /*
     * здесь создаётся маленький граф для примеров трансформации
     */

    const std::string NUCLS = "ACGT";
    std::srand(std::time(nullptr));
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,3);

    std::ostringstream oss;
    for (size_t j = 0; j < genome_length; ++j) {
        oss << NUCLS[distribution(generator)];
    }
    std::string genome = oss.str();
    //std::cout << genome << "\n";
    unsigned k = (unsigned)g.k();

    std::vector<std::string> reads;

    for (size_t i = 0; i < 1/*3500*/; ++i) {
        size_t start = rand() % (genome.length() - k - 1000);
        size_t len = 60;//rand() % (genome.length() - k - 100) + (k + 100);
        //for (int j = 0; j < 100; ++j)
            reads.push_back(genome.substr(start, len));
    }

    INFO("reads generated");
    // reads != kmers

    // вообще-то list не нужен, но я не знаю, как иначе измерить cardinality
    io::ReadStreamList<io::SingleReadSeq> read_streams(std::move(StreamWrapper(reads)));

    INFO("read stream is created");
    // ReadStream<ReadType>. из них можно читать риды

    // io::SingleReadSeq common/io/reads/single_read.hpp
    // хранит последовательность нуклеитидов какого-то рида +
    // смещение влево и вправо от исходной последовательности

/*
 config::debruijn_config::construction params:
   arly_tip_clipper early_tc;
   bool keep_perfect_loops;
   unsigned read_cov_threshold;
   size_t read_buffer_size;
 */

    kmers::DeBruijnExtensionIndex<> ext_index(k); // common/kmer_index/extension_index/kmer_extension_index.hpp

    // по умолчанию ext_index::storing_type = InvertableStoring - common/kmer_index/ph_map/storing_traits.hpp
   using storing_type = decltype(ext_index)::storing_type;

// 1. CoverageFilter

    unsigned k_plus_one = ext_index.k() + 1;
    unsigned read_coverage_threshold = 1; // а чему оно должно быть равно в норме?

    using KmerFilter = kmers::StoringTypeFilter<storing_type>;


    rolling_hash::SymmetricCyclicHash<rolling_hash::NDNASeqHash> hasher(k_plus_one);

    size_t kmers_cardinality = EstimateCardinalityUpperBound(k_plus_one, read_streams, hasher, KmerFilter()); //wtf почему опять нужен stream.size()??
    std::cout << "kmers_cardinality " << kmers_cardinality << "\n";
    kmers_cardinality = 10;

    INFO("Estimated k-mers cardinality");

    qf::cqf cqf(kmers_cardinality); // common/adt/cqf.hpp

    INFO("Building k-mer coverage histogram");
    FillCoverageHistogram(cqf, k_plus_one, hasher, read_streams, read_coverage_threshold, KmerFilter());

    // вот до сюда дошли без ошибок

    // Replace input streams with wrapped ones
    read_streams = io::CovFilteringWrap(std::move(read_streams), k_plus_one, hasher, cqf, read_coverage_threshold);
// wtf я добавила #include "paired_reads" в coverage_filtering_read_wrapper.hpp, и оно заработало, но spades прекрасно работает без этого???
    /*
    // io::CovFilteringWrap common/io/reads/coverage_filtering_read_wrapper.hpp -
    // обёртка для ReadStream / ReadStreamList, в которых по какому-то там признаку отфильтрованы риды

*/
    INFO("coverage filtering finished");

// 2. KMerCounting - тут есть Splitter, он режет риды на камеры
    unsigned n_threads = (unsigned)read_streams.size();

    using Splitter =  kmers::DeBruijnReadKMerSplitter<io::SingleReadSeq, kmers::StoringTypeFilter<storing_type>>;

    auto work_dir = fs::impl::make_temp_dir("", "work_dir");
    size_t read_buffer_size = 0;   // по умолчанию 0 в конструкторе Splitter
    kmers::KMerDiskCounter<RtSeq> counter(work_dir,
                                          Splitter(work_dir, k_plus_one, read_streams, read_buffer_size));

    std::cout << "k+1 " << k_plus_one << "\n";
    kmers::KMerDiskStorage<RtSeq> kmers = counter.Count(10 * n_threads, n_threads);

    // kmers::KMerDiskCounter - common/kmer_index/kmer_mph/kmer_index_builder.hpp
    // kmers::DeBruijnReadKMerSplitter - common/kmer_index/kmer_mph/kmer_splitters.hpp
    // kmers::StoringTypeFilter - common/kmer_index/ph_map/storing_traits.hpp
    // kmers::KMerDiskStorage<RtSeq> - common/kmer_index/kmer_mph/kmer_index_builder.hpp. есть конструктор перемещения и присваивания через перемещение



// 3. ExtensionIndexBuilder
// тут тоже по умолчанию read_buffer_size = 0
    kmers::DeBruijnExtensionIndexBuilder().BuildExtensionIndexFromKPOMers(work_dir,
                                                                          ext_index,
                                                                          kmers,
                                                                          read_streams.size(),
                                                                          read_buffer_size);
    INFO("extension index is built");

    std::cout << ext_index.size() << "\n";
// 4. GraphCondenser
    debruijn_graph::DeBruijnGraphExtentionConstructor<debruijn_graph::Graph> graph_constructor(g,ext_index); // common/assembly_graph/construction/debruijn_graph_constructor.hpp
    INFO("constructor created");

    bool keep_perfect_loops = 0; // что именно это значит?
    graph_constructor.ConstructGraph(keep_perfect_loops);
    INFO("graph constructed");

// 5. PHMCoverageFiller - надо?

// ?. Сохранение графа
    spades_example::SaveToGFA(g, "random_graph.gfa");
    INFO("graph is saved");


}
} // spades_example


void CreateConsoleLogger(const std::filesystem::path& log_fn="") {
    using namespace logging;
    logger *lg = create_logger(exists(log_fn) ? log_fn : "");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

int main() {
    CreateConsoleLogger();
    debruijn_graph::Graph g(55);
    // genome_length should be bigger (much bigger actually) then k
    spades_example::CreateRandomGraph(g, 10000);
}