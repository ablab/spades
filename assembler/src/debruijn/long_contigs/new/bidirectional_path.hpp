/*
 * bidirectional_path.h
 *
 *  Created on: Nov 14, 2011
 *      Author: andrey
 */

#ifndef BIDIRECTIONAL_PATH_H_
#define BIDIRECTIONAL_PATH_H_

#include "../new_debruijn.hpp"

using debruijn_graph::Graph;
using debruijn_graph::EdgeId;
using debruijn_graph::VerexId;


class BidirectionalPath;


class LoopDetectorData {

typedef std::map <EdgeId, double> AltenativeMap;

protected:
    size_t iteration;

    AltenativeMap alternatives;

public:
    LoopDetectorData(size_t i): iteration(i), alternatives()  {
    }

    LoopDetectorData(): alternatives()  {
    }

    LoopDetectorData(const LoopDetectorData& d) {
        iteration = d.iteration;
        alternatives.insert(d.alternatives.begin(), d.alternatives.end());
    }

    size_t getIteration() const {
        return iteration;
    }

    const AltenativeMap& getAlternatives() const {
        return alternatives;
    }

    void addAlternative(EdgeId e, double w = 1) {
        alternatives.insert(std::make_pair(e, w));
    }

    void clear() {
        alternatives.clear();
        iteration = 0;
    }

    bool operator==(const LoopDetectorData& d) const {
        if (alternatives.size() != d.alternatives.size()) {
            return false;
        }

        auto iter2 = d.alternatives.begin();
        for (auto iter1 = alternatives.begin(); iter2 != d.alternatives.end() && iter1 != alternatives.end(); ++iter1, ++iter2) {
            if (iter1->first != iter2->first || iter1->second != iter2->second) {
                return false;
            }
        }

        return true;
    }
};


class LoopDetector {

protected:
    Graph& g;
    BidirectionalPath& parent;

    size_t currentIteration;
    size_t loopLimit;

    LoopDetectorData * current;
    std::multimap <EdgeId, LoopDetectorData* > data;

public:
    LoopDetector(Graph& g_, BidirectionalPath& p_, size_t loopLimit_);

    void addAlternative(EdgeId e, double w = 1);

    void selectEdge(EdgeId e, double weight = 1);

    void clear();

    ~LoopDetector();

    size_t loopEdges();

    size_t loopLength();

    bool pathIsLoop() const;

    size_t lastLoopCount() const;

    bool isCycled() const;

    size_t edgesToRemove(bool fullRemoval = false);

    void removeLoop(size_t loopCount = 0, bool fullRemoval = false);

    bool LoopBecameStable() const;

    size_t CountLoopExits(BidirectionalPath& path, EdgeId e, LoopDetector& detector, bool forward);

    EdgeId FindFirstFork(BidirectionalPath& path, EdgeId e, LoopDetector& detector, bool forward);

    EdgeId GetForwardFork(Graph& g, EdgeId e);

    EdgeId GetBackwardFork(Graph& g, EdgeId e);

    bool EdgesMakeShortLoop(Graph& g, EdgeId e1, EdgeId e2);

    EdgeId IsEdgeInShortLoopForward(Graph& g, EdgeId e);

    EdgeId IsEdgeInShortLoopBackward(Graph& g, EdgeId e);

    size_t GetMaxExitIteration(EdgeId loopEdge, EdgeId loopExit, LoopDetector& detector, std::pair<size_t, size_t> iterRange);

    size_t GetFirstExitIteration(EdgeId loopEdge, EdgeId loopExit, LoopDetector& detector,  std::pair<size_t, size_t> iterRange, double coeff);

    void print() const {
        INFO("== Detector data ==");
        for (auto iter = data.begin(); iter != data.end(); ++iter) {
            INFO("Edge " << g.length(iter->first));

            const LoopDetectorData::AltenativeMap& alts = iter->second->getAlternatives();
            for(auto alt = alts.begin(); alt != alts.end(); ++alt) {
                INFO("Edge " << g.length(alt->first) << ", weight " << alt->second);
            }
        }
    }

};


class BidirectionalPath {

protected:
	Graph& g;

	//Unique ID
	size_t id;

	// Edges: e1 e2 ... eN
	std::deque <EdgeId> data;

	// Length from beginning of i-th edge to path end for forward directed path:
	// L(e2 + ... + eN) ... L(eN), 0
	// Length from beginning of the path to the end of i-th edge
	// L(e1), L(e1 + e2) ... L(e1 + ... + eN)
	std::deque <size_t> cumulativeLength;

	// e1 - gap1 - e2 - ... - gap(N-1) - eN
	std::deque <size_t> gapLength;

	// L(e1 + ... + eN)
	size_t totalLength;

	// Defines meaning of cumulative length
	bool direction;

	// Cycle analyzer
	LoopDetector loopDetector;

protected:
    void recountLengths() {

        cumulativeLength.clear();
        size_t currentLength = 0;

        if (direction) {
            for(auto iter = data.rbegin(); iter != data.rend(); ++iter) {
                currentLength += g.length(*iter);
                cumulativeLength.push_front(currentLength);
            }
        } else {
            for(auto iter = data.begin(); iter != data.end(); ++iter) {
                cumulativeLength.push_back(currentLength);
                currentLength += g.length(*iter);
            }
        }

        totalLength = currentLength;
    }

    void increaseLengths(size_t length, bool direct = direction) {
        for(auto iter = cumulativeLength.begin(); iter != cumulativeLength.end(); ++iter) {
            *iter += length;
        }

        if (direct) {
            cumulativeLength.push_back(length);
        } else {
            cumulativeLength.push_front(0);
        }

        totalLength += length;
    }

    void decreaseLengths(bool direct = direction) {
        size_t length = direct ? g.length(data.back()) : g.length(data.front());
        for(auto iter = cumulativeLength.begin(); iter != cumulativeLength.end(); ++iter) {
            *iter -= length;
        }

        if (direct) {
            cumulativeLength.pop_back();
        } else {
            cumulativeLength.pop_front();
        }

        totalLength -= length;
    }

public:
	BidirectionalPath(Graph g_): g(g_), data(), cumulativeLength(), gapLength(), direction(true), loopDetector(g_, *this) {
	}

	size_t size() const {
	    return data.size();
	}

	size_t length() const {
	    return totalLength;
	}

	//Access methods
	EdgeId operator[](size_t index) const {
	    return data[index];
	}

	EdgeId at(size_t index) const {
	    return data[index];
	}

	size_t lengthAt(size_t index) const {
	    return cumulativeLength[index];
	}

	size_t gapAt(size_t index) const {
	    return gapLength[index];
	}

	bool getDirection() const {
	    return direction;
	}

	size_t getId() {
	    return id;
	}

	EdgeId head() const {
	    return direction ? data.back() : data.front();
	}

	EdgeId back() const {
	    return data.back();
	}

	EdgeId front() const {
	    return data.front();
	}

	LoopDetector& getLoopDetector() {
	    return loopDetector;
	}

	//Modification methods
	void setDirection(bool d) {
	    direction = d;
	}

	void pushBack(EdgeId e, size_t gap = 0) {
	    data.push_back(e);
	    gapLength.push_back(gap);
	    increaseLengths(g.length(e), true);
	}

	void pushFront(EdgeId e, size_t gap = 0) {
	    data.push_front(e);
	    gapLength.push_front(gap);
	    increaseLengths(g.length(e), false);
	}

	void popBack() {
	    data.pop_back();
	    gapLength.pop_back();
	    decreaseLengths(true);
	}

	void popFront() {
	    data.pop_front();
	    gapLength.pop_front();
	    decreaseLengths(false);
	}


	void push(EdgeId e, size_t gap = 0) {
	    if (direction) {
	        pushBack(e, gap);
	    } else {
	        pushFront(e, gap);
	    }
	}

	bool pop() {
	    if (data.empty()) {
	        return false;
	    }
	    if (direction) {
	        popBack();
	    } else {
	        popFront();
	    }
	    return true;
	}

	void setId(size_t uid) {
	    id = uid;
	}

    BidirectionalPath(Graph g_, EdgeId startingEdge) {
        BidirectionalPath(g_);
        push(startingEdge);
        id = g.int_id(startingEdge);
    }

};


class PathContainer {
private:
	std::vector <BidirectionalPath*> data;

public:
	PathContainer() {
		data.clear();
	}

    BidirectionalPath* get(size_t index) const {
        return data[index];
    }

    BidirectionalPath* getConjugate(size_t index) const {
        return data[(index & 1 == 0) ? index + 1 : index - 1];
    }

    void reserve(size_t size) {
        data.reserve(size);
    }

    bool addPair(BidirectionalPath* p, BidirectionalPath* cp) {
        if (p->size() != cp->size() || p->length() != cp->length()) {
            return false;
        }

        data.push_back(p);
        data.push_back(cp);
        return true;
    }

    bool addPair(BidirectionalPath& p, BidirectionalPath& cp) {
        if (p.size() != cp.size() || p.length() != cp.length()) {
            return false;
        }

        BidirectionalPath * np = new BidirectionalPath(p);
        BidirectionalPath * ncp = new BidirectionalPath(cp);

        data.push_back(np);
        data.push_back(ncp);
        return true;
    }
};


LoopDetector::LoopDetector(Graph& g_, BidirectionalPath& p_, size_t loopLimit_): g(g_), parent(p_), loopLimit(loopLimit_), currentIteration(0) {
  current = new LoopDetectorData(currentIteration);
}

void LoopDetector::addAlternative(EdgeId e, double w = 1) {
  current->addAlternative(e, w);
}


void selectEdge(EdgeId e, double weight = 1) {
  data.insert(std::make_pair(e, current));
  current = new LoopDetectorData(++currentIteration);
}


void clear() {
  for (auto iter = data.begin(); iter != data.end(); ++iter) {
      delete iter->second;
  }

  data.clear();
  current->clear();
}


~LoopDetector() {
  clear();
  delete current;
}


size_t loopEdges() const {
  EdgeId e = parent.head();

  if (data.count(e) <= 1) {
      return 0;
  }

  auto iter = data.upper_bound(e);
  --iter;
  size_t loopSize = iter->second->getIteration();
  --iter;
  loopSize -= iter->second->getIteration();

  return loopSize;
}


size_t loopLength() const {
  size_t edges = loopEdges();
  size_t length = 0;

  if (parent.getDirection()) {
      for (int i = parent.size() - edges; i < parent.size(); ++i) {
          length += g.length(parent[i]);
      }
  } else {
      for (int i = 0; i < edges; ++i) {
          length += g.length(parent[i]);
      }
  }
  return length;
}


bool pathIsLoop() const {
  size_t edges = loopEdges();

  if (edges <= parent.size()) {
      return false;
  }

  for (int i = 1; i <= edges; ++i) {
      EdgeId e = parent[parent.size() - i];
      for (int j = parent.size() - i - edges; j >= 0; j -= edges) {
          if (parent[j] != e) {
              return false;
          }
      }
  }
  return true;
}


size_t lastLoopCount() const {
  size_t edges = loopEdges();
  EdgeId e = parent.head();
  size_t count = 0;

  int i = parent.getDirection() ? parent.size() - 1 - edges : 0;
  int delta = parent.getDirection() ? -edges : edges;

  while (parent[i] == e) {
      ++count;
      i += delta;
  }

  return count;
}


bool isCycled() const {
  return lastLoopCount() >= loopLimit;
}


size_t edgesToRemove(bool fullRemoval = false) {
  size_t edges = loopEdges();
  size_t count = lastLoopCount();
  bool onlyCycle = pathIsLoop();

  if (onlyCycle || parent.size() <= count * edges + 1) {
      return  parent.size() - edges;
  }

  if (fullRemoval) {
      return count * edges;
  } else {
      return (count - 1) * edges + 1;
  }
}

void removeLoop(size_t loopCount = 0, bool fullRemoval = false) {
  if (loopCount == 0) {
      loopCount = lastLoopCount();
  }

  size_t edgesToRemove = edgesToRemove(fullRemoval);

  for(size_t i = 0; i < edgesToRemove; ++i) {
      parent.pop();
  }
}

bool LoopBecameStable() const {
  Edge e = parent.head();

  if (detector.data.count(e) < 2) {
      DETAILED_INFO("Loop still unstable");
      return false;
  }
  auto iter = detector.data.upper_bound(e);
  auto last = --iter;
  auto prev = --iter;

  bool res = prev->second == last->second;

  if (res) {
      DETAILED_INFO("Loop became stable");
  } else {
      DETAILED_INFO("Loop still unstable");
  }
  return res;
}

size_t CountLoopExits(BidirectionalPath& path, EdgeId e, LoopDetector& detector, bool forward) {
  size_t loopSize = CountLoopEdges(e, detector);
  size_t exits = 0;
  int start = forward ? path.size() - 1 : loopSize - 1;
  int end = forward ? path.size() - loopSize : 0;

  for (int i = start; i >= end; --i) {
      LoopDetectorData& data = detector.data.find(path[i])->second;

      exits += data.weights.size() - 1;
  }
  return exits;
}

EdgeId FindFirstFork(BidirectionalPath& path, EdgeId e, LoopDetector& detector, bool forward) {
  size_t loopSize = CountLoopEdges(e, detector);
  int start = forward ? path.size() - 1 : loopSize - 1;
  int end = forward ? path.size() - loopSize : 0;

  for (int i = start; i >= end; --i) {
      LoopDetectorData& data = detector.data.find(path[i])->second;

      if (data.weights.size() == 2) {
          return path[i];
      }
  }
  return 0;
}

EdgeId GetForwardFork(Graph& g, EdgeId e) {
  VertexId v = g.EdgeStart(e);
  if (g.OutgoingEdgeCount(v) != 2) {
      return 0;
  }
  std::vector<EdgeId> edges = g.OutgoingEdges(v);
  if (edges[1] == e) {
      return edges[0];
  } else {
      return edges[1];
  }
}

EdgeId GetBackwardFork(Graph& g, EdgeId e) {
  VertexId v = g.EdgeEnd(e);
  if (g.IncomingEdgeCount(v) != 2) {
      return 0;
  }
  std::vector<EdgeId> edges = g.IncomingEdges(v);
  if (edges[1] == e) {
      return edges[0];
  } else {
      return edges[1];
  }
}

bool EdgesMakeShortLoop(Graph& g, EdgeId e1, EdgeId e2) {
  return g.EdgeStart(e1) == g.EdgeEnd(e2) && g.EdgeStart(e2) == g.EdgeEnd(e1);
}

EdgeId IsEdgeInShortLoopForward(Graph& g, EdgeId e) {
  VertexId v = g.EdgeEnd(e);
  auto edges = g.OutgoingEdges(v);
  EdgeId result = 0;

  for (auto edge = edges.begin(); edge != edges.end(); ++edge) {
      if (g.EdgeEnd(*edge) == g.EdgeStart(e)) {
          result = *edge;
      }
  }

  if (g.OutgoingEdgeCount(v) == 1 && result != 0) {
      INFO("Seems no fork backward: edge " << g.length(e) << ", loops with " << g.length(result) << ". " << g.OutgoingEdgeCount(v));
  }

  return result;
}

EdgeId IsEdgeInShortLoopBackward(Graph& g, EdgeId e) {
  VertexId v = g.EdgeStart(e);
  auto edges = g.IncomingEdges(v);
  EdgeId result = 0;

  for (auto edge = edges.begin(); edge != edges.end(); ++edge) {
      if (g.EdgeStart(*edge) == g.EdgeEnd(e)) {
          result = *edge;
      }
  }

  if (g.IncomingEdgeCount(v) == 1 && result != 0) {
      INFO("Seems no fork backward: edge " << g.length(e) << ", loops with " << g.length(result) << ". " << g.IncomingEdgeCount(v));
  }

  return result;
}

size_t GetMaxExitIteration(EdgeId loopEdge, EdgeId loopExit, LoopDetector& detector, std::pair<size_t, size_t> iterRange) {
  auto range = detector.data.equal_range(loopEdge);

  size_t maxIter = 0;
  double maxWeight = 0;
  for (auto iter = range.first; iter != range.second; ++iter) {
      double w = iter->second.weights[loopExit];
      if (w > maxWeight &&
              iter->second.iteration >= iterRange.first && iter->second.iteration <= iterRange.second) {

          maxIter = iter->second.iteration;
          maxWeight = w;
      }
  }
  return maxIter;
}

size_t GetFirstExitIteration(EdgeId loopEdge, EdgeId loopExit, LoopDetector& detector,
      std::pair<size_t, size_t> iterRange, double coeff = lc_cfg::get().es.priority_coeff) {
  auto range = detector.data.equal_range(loopEdge);

  size_t maxIter = std::numeric_limits<size_t>::max();
  for (auto iter = range.first; iter != range.second; ++iter) {
      if (iter->second.weights[loopExit] * coeff > iter->second.weights[loopEdge] && maxIter > iter->second.iteration &&
              iter->second.iteration >= iterRange.first && iter->second.iteration <= iterRange.second) {

          maxIter = iter->second.iteration;
      }
  }
  return maxIter;
}






#endif /* BIDIRECTIONAL_PATH_H_ */
