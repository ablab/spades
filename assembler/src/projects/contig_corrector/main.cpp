#include "sequence_corrector/sequence_corrector.hpp"
#include "contig_replacer/contig_replacer.hpp"
#include "error_analyzer/error_analyzer.hpp"
#include "sequence_clusterer/sequence_clusterer.hpp"

#include "common/toolchain/utils.hpp"
#include "utils/segfault_handler.hpp"

#include <clipp/clipp.h>

#include <iostream>

int main(int argc, char *argv[]) {
    utils::segfault_handler sh;
    toolchain::create_console_logger(logging::L_INFO);
    srand(42);
    srandom(42);

    try {
        enum class Action {
            replace,
            correct,
            analyze,
            clustering,
        } action;

        auto cli = (
            (clipp::command("replace").set(action, Action::replace), contig_replacer::GetCLI()) |
            (clipp::command("correct").set(action, Action::correct), sequence_corrector::GetCLI()) |
            (clipp::command("analyze").set(action, Action::analyze), error_analyzer::GetCLI()) |
            (clipp::command("clustering").set(action, Action::clustering), sequence_clusterer::GetCLI())
        );

        auto result = clipp::parse(argc, argv, cli);
        if (!result) {
            if (result.begin() != result.end()) {
                std::cerr << "Args parsed as: <arg, index>\n";
                for (auto const & res : result)
                    std::cerr << res.arg() << ' ' << res.index() << ' ' << '\n';
            }
            std::cerr << clipp::make_man_page(cli, argv[0]) << std::endl;
            return -1;
        }

        switch (action) {
        case Action::replace: 
            return contig_replacer::main();
        case Action::correct: 
            return sequence_corrector::main();
        case Action::analyze: 
            return error_analyzer::main();
        case Action::clustering: 
            return sequence_clusterer::main();
        default: assert(false && "unreachable");
        }
    } catch (const std::string &s) {
        std::cerr << s << std::endl;
        return EINTR;
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return EINTR;
    } catch (...) {
        std::cerr << "unknown exception caught" << std::endl;
        return EINTR;
    }
    return 0;
}
