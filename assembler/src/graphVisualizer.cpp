#import "graphVisualizer.hpp"

namespace gvis {

void startGraphRecord(string name) {
	cout << "digraph " << name << " {" << endl;
}

void endGraphRecord() {
	cout << "}" << endl;
}

}
