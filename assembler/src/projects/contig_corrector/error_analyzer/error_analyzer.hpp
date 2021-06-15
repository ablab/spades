#pragma once

#include <clipp/clipp.h>

namespace error_analyzer {

int main(int argc, char * argv[]);

clipp::group GetCLI();

} // namespace error_analyzer
