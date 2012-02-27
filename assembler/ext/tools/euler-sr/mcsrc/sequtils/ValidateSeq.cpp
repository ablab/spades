#include "SeqReader.h"
#include "DNASequence.h"
#include <string>
using namespace std;

int main(int argc, char* argv[]) {

	string inName = argv[1];
	ifstream in;
	openck(inName, in, std::ios::in);
	string title;
	while(in) {
		if (in and in.peek() == '>') {
			string title;
			std::getline(in, title);
		}
		else {
			char c;
			while (in and ((c = in.get()) != '>')) {
				if ( (c < 'A' or c > 'Z') and
						 (c < 'a' or c > 'z') and
						 c != '\n' and c != '\t' and c != ' ') {
					cout << "read: " << title << " is bad." << endl;
				}
			}
		}
	}

}
