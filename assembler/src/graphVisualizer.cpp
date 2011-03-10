#include "graphVisualizer.hpp"

namespace gvis {

void startGraphRecord(ostream &out, string &name) {
	out << "digraph " << name << " {" << endl;
}

void endGraphRecord(ostream &out) {
	out << "}" << endl;
}

}
