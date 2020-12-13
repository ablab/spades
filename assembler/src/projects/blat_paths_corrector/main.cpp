#include "sequence_corrector.hpp"

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
        } action;

        auto cli = (
            clipp::command("replace").set(action, Action::replace) | 
            (
                clipp::command("correct").set(action, Action::correct),
                sequence_corrector::GetCLI()
            )
        );

        if (!clipp::parse(argc, argv, cli)) {
            std::cout << clipp::make_man_page(cli, argv[0]) << '\n';
            return -1;
        }

        switch (action) {
        case Action::replace: 
            return -1;
        case Action::correct: 
            return sequence_corrector::main(argc, argv);
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
