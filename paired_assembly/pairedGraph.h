#ifndef PAIREDGRAPH_H_
#define PAIREDGRAPH_H_
#include <vector>
typedef char *sequence;

struct Node {
    sequence upperString;
    sequence lowerString;
    short upperSize;
    short lowerSize;
    Node* neighbours;
//may be only d or only delta
    char delta;
    char d;
//bounds for d?
};
//TEMPORARY
typedef std::vector<Node *> graph;
typedef std::vector<Node *> nodelist;
typedef long long Kmer;
typedef sequence Read;

int isDeleted(Node* A);
int deleteNode(Node* A);

int addNode(Node* A);
int addEdgesFromNode(Node* A);

int addEdge(Node* A , Node * B);

//unused? 
int removeEdge(Node* A, Node* B);

int appendNode(Node* A, Node* B);
//Добавляет B к A.
//Храним ли мы логи?
//возвращаем удалось ли смерджить
//интеллектуально выделять память.

int merge(Node * A, Node * B);

nodelist findKmer(Kmer k);

nodelist findRead(Read r);

//В частности, удаляем "удаленные" вершины; возможно она не graph а void.
graph rearrange();

#endif /*PAIREDGRAPH_H_*/
