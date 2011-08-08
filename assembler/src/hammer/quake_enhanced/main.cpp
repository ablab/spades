#include <iostream>
#include "common/logging.hpp"
#include "hammer/quake_enhanced/options.hpp"
DECL_PROJECT_LOGGER("q")

namespace quake_enhanced {

DECL_LOGGER("main")

}

int main(int argc, char **argv) {
  quake_enhanced::Options opts(argc, argv);
  if (!opts.valid) {
    std::cout << opts.help_message;
    return 1;
  }
  return 0;
}
