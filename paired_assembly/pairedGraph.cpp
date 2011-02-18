#include "pairedGraph.h"

int isDeleted(Node * A) {
	return A->upperSize == 0;
}
int deleteNode(Node * A) {
	return A->upperSize = 0;
}