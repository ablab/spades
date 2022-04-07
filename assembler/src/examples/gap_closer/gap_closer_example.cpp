#include "projects/spades/gap_closer.hpp"
#include "common/assembly_graph/core/graph.hpp"
#include "common/utils/logger/log_writers.hpp"
#include "common/pipeline/graph_pack.hpp"

#include <filesystem>
#include <string>
#include <unordered_set>

#include <gtest/gtest.h>

void AssertEdges(const debruijn_graph::Graph& g1, const debruijn_graph::Graph& g2) {
    INFO("Asserting edges");
    std::unordered_set<std::string> edges1;
    //constiterator, и вообще посмотреть разные итераторы
    for (auto it = g1.SmartEdgeBegin(); !it.IsEnd(); ++it) {
        edges1.insert(g1.EdgeNucls(*it).str());
    }
    std::unordered_set<std::string> edges2;
    for (auto it = g2.SmartEdgeBegin(); !it.IsEnd(); ++it) {
        edges2.insert(g2.EdgeNucls(*it).str());
    }

    ASSERT_EQ(edges1.size(), edges2.size());
    for (auto it = edges1.begin(); it != edges1.end(); ++it) {
        ASSERT_TRUE(edges2.count(*it) > 0);
    }
}

void load_and_save_graph_pack() {
    //пути локальные убрать, чтобы универсально было
    //инструкцию писать из расчета на то, что я хотела бы прочитать, когда начинала это всё
    std::filesystem::path load_from = "/Bmo/aberdichevskaya/algorithmic-biology/assembler/tmp/K55/saves";
    std::filesystem::path work_dir = "/Bmo/aberdichevskaya/algorithmic-biology/assembler/tmp_dir/tmp";
    std::filesystem::path binary_reads_dir = "/Bmo/aberdichevskaya/algorithmic-biology/assembler/tmp/.bin_reads";

    std::filesystem::path output_dir = "/Bmo/aberdichevskaya/algorithmic-biology/assembler/src/examples/gap_closer/tmp";

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

void create_console_logger(std::filesystem::path const& log_fn) {
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