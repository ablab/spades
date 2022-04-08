#include "projects/spades/gap_closer.hpp"
#include "common/assembly_graph/core/graph.hpp"
#include "common/utils/logger/log_writers.hpp"
#include "common/pipeline/graph_pack.hpp"

#include <filesystem>
#include <string>
#include <unordered_set>

#include <gtest/gtest.h>

void AssertEdges(const debruijn_graph::Graph& g1, const debruijn_graph::Graph& g2) {
    /*
     * Edges of de Bruijn graph present k-1 mers in a form of std::string. Thus, edges of the graph
     * can be compared as a set of std::string.
     *
     * There are several iterators to go through edges in debruijn_graph::Graph. They are declared in
     * common/assembly_graph/core/observable_graph.hpp file.
     */

    INFO("Asserting edges");
    std::unordered_set<std::string> edges1;
    for (auto it = g1.ConstEdgeBegin(); !it.IsEnd(); ++it) {
        edges1.insert(g1.EdgeNucls(*it).str());
    }
    std::unordered_set<std::string> edges2;
    for (auto it = g2.ConstEdgeBegin(); !it.IsEnd(); ++it) {
        edges2.insert(g2.EdgeNucls(*it).str());
    }

    ASSERT_EQ(edges1.size(), edges2.size());
    for (auto it = edges1.begin(); it != edges1.end(); ++it) {
        ASSERT_TRUE(edges2.count(*it) > 0);
    }
}

void load_and_save_graph_pack() {
    /*
     * This example shows how to use one of the SPAdes stages (in this case it is Gap Closer stage)
     * separately from the SPAdes pipeline.
     *
     *
     */
    //пути локальные убрать, чтобы универсально было
    //инструкцию писать из расчета на то, что я хотела бы прочитать, когда начинала это всё
    std::filesystem::path load_from = "/Bmo/aberdichevskaya/algorithmic-biology/assembler/tmp/K55/saves";
    std::filesystem::path work_dir = "/Bmo/aberdichevskaya/algorithmic-biology/assembler/tmp_dir/tmp";
    std::filesystem::path binary_reads_dir = "/Bmo/aberdichevskaya/algorithmic-biology/assembler/tmp/.bin_reads";

    std::filesystem::path output_dir = "/Bmo/aberdichevskaya/algorithmic-biology/assembler/src/examples/gap_closer/tmp";

    /*
     * GraphPack is a data structure, which stores information about assembly de Bruijn graph. It is declared
     * in common/pipeline/graph_pack.hpp file.
     *
     * Data storage is organized in a such way, that only one instance of GraphPack may exist during the runtime.
     *
     * Access to information about the graph is carried out via get<T>() and get_mutable<T>() functions,
     * where T is a name of data structure in which this information is stored.
     *
     * GraphPack has 3 obligatory arguments: size_t k - the length of reads, std::filesystem::path work_dir -
     * some directory, which GraphPack uses during the runtime, and size_t lib_count (which is what?).
     */
    graph_pack::GraphPack gp(55, work_dir, 1); //std::vector<std::string>(0), 55, 0, 0, true);
    INFO("created graph pack");

    debruijn_graph::GapClosing gapClosingStage("gap_closer_example");
    INFO("created gap closing stage");

    create_directory(binary_reads_dir);

    gapClosingStage.load(gp, load_from, "ec_threshold_finder");
    //вывести какую-то инфу (размер графа, зааттачены ли индексы, что индекс какой-нибудь не пустой)

    INFO("graph is loaded");

    gapClosingStage.run(gp, "");

    INFO("gap closing finished");

    remove_all(output_dir);
    create_directory(output_dir);

    gapClosingStage.save(gp, output_dir, "gap_closer");

    INFO("graph is saved");

    const auto& g1 = gp.get<debruijn_graph::Graph>();


    gapClosingStage.load(gp, load_from, "early_gapcloser");
    const auto& g2 = gp.get<debruijn_graph::Graph>();

    AssertEdges(g1, g2);

    INFO("edges of graphs are equal");
}

void load_and_save_graph() {
    //загрузить и выгрузить граф (.grseq) проверить, что норм это происходит

    //можно загружать/выгружать сам граф, а не графпак

    //сделать ещё загрузку не через stage, а просто граф
}

void create_console_logger(std::filesystem::path const& log_fn = "") {
/*
 * To track the process of program execution it is useful to create a console logger.
 * In SPAdes logger is declared in common/utils/logger/logger.hpp file and logger writer is declared
 * in common/utils/logger/log_writters.hpp. To use logger in the project common/utils/logger/log_writters.hpp
 * should be included (in this case common/utils/logger/logger.hpp is including transitive) and link
 * utils library in CMake file.
 *
 * Logs will be written in log_fn file. (а собственно где этот файл?). In case if log_fn is an empty path,
 * (а что собственно случится тогда?).
 *
 * After logger creation it is convenient to use TRACE, DEBUG, INFO, WARN and ERROR macro.
 */
    using namespace logging;
    logger *lg = create_logger(exists(log_fn) ? log_fn : "");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

int main() {
    create_console_logger("log.txt");
    INFO("Start of the example");

    ///
    std::vector<std::filesystem::path> cfg_fns = {"/Bmo/aberdichevskaya/algorithmic-biology/assembler/tmp/K55/configs/config.info",
                                        "/Bmo/aberdichevskaya/algorithmic-biology/assembler/tmp/K55/configs/mda_mode.info"};
    cfg::create_instance(cfg_fns);
    ////

    load_and_save_graph_pack();

    return 0;
}