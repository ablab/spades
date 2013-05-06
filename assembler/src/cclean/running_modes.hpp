#ifndef RUNNING_MODES_H_
#define RUNNING_MODES_H_

#include <ostream>
#include "database.hpp"
#include "io/ireadstream.hpp"

void exactMatch(std::ostream& output, std::ostream& bed, ireadstream * input, const Database * data);
void alignment(std::ostream& output, std::ostream& bed, ireadstream * input, const Database * data);
void exactAndAlign(std::ostream& output, std::ostream& bed, ireadstream * input, const Database * data);

#endif /* RUNNING_MODES_H_ */
