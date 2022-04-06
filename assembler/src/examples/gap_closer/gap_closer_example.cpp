#include "projects/spades/gap_closer.hpp"
#include "common/utils/logger/log_writers.hpp"
#include "common/pipeline/graph_pack.hpp"

#include <filesystem>

void load_and_save_graph_pack() {
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

    //compare graphs
    //AssertGraph (test/debruijn/test_utils.cpp : 188)
}

void load_and_save_graph() {
    //загрузить и выгрузить граф (.grseq)  проверить, что норм это происходит

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