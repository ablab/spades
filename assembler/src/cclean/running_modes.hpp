#ifndef RUNNING_MODES_H_
#define RUNNING_MODES_H_

#include <ostream>
#include <unordered_map>
#include <string>
#include "additional.cpp"

class ireadstream;
namespace cclean {
class AdapterIndex;
}

void ExactAndAlign(std::ostream &aligned_output,
                   std::ostream &bed, ireadstream *input, std::ostream &output,
                   const std::string &db, const cclean::AdapterIndex &index,
                   const additional::WorkModeType &mode,
                   const std::unordered_map<std::string, std::string> &options);

#endif /* RUNNING_MODES_H_ */
