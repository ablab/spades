#include "config_struct.hpp"

#include "openmp_wrapper.h"

#include <yaml-cpp/yaml.h>
#include <string>

namespace corrector {
void load(corrector_config& cfg, const std::string &filename) {
    YAML::Node config = YAML::LoadFile(filename);
    cfg.dataset.load(config["dataset"].as<std::string>());
    // WTF: Why extra spaces? Fix the coding style everywhere!
    // Re: Because of official google-SPAdes style (style file for eclipse by AntonB, ext/eclipse/gsgc.xml)
    // Re: I fixed extra spaces, but if this file is not real SPAdes code style description then where it is?
    // WTF: The only style guide is located at http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml
    // Everything else is an approximation. Period. Even more, in this file you can easily find:
    // <setting id="org.eclipse.cdt.core.formatter.insert_space_after_opening_angle_bracket_in_template_arguments" value="do not insert"/>
    // So, maybe you're not using it after all?
    cfg.work_dir = config["work_dir"].as<std::string>(".");
    cfg.output_dir = config["output_dir"].as<std::string>(".");
    cfg.max_nthreads = config["max_nthreads"].as<unsigned>();
    // Fix number of threads according to OMP capabilities.
    cfg.max_nthreads = std::min(cfg.max_nthreads, (unsigned) omp_get_max_threads());
    // WTF: Add proper serialization / deserialization of the enum. Learn about enum classes. There are plenty of examples around.
    std::string strategy_str = config["strategy"].as<std::string>(".");
    if (strategy_str == "all_reads")
        cfg.strat = strategy::all_reads;
    else if (strategy_str == "majority_only")
        cfg.strat = strategy::majority_only;
    else if (strategy_str == "not_started")
        cfg.strat= strategy::not_started;
    else
        //default strat
        cfg.strat= strategy::mapped_squared;
    cfg.bwa = config["bwa"].as<std::string>(".");
// Inform OpenMP runtime about this :)
    omp_set_num_threads(cfg.max_nthreads);
}
}
