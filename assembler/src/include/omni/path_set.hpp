#pragma once
#include <cmath>
#include <map>
#include <limits>
#include <xmath.h>

namespace omnigraph{

/**
 * PathSet: the basic data structure for the path-set graph
 **/
template<typename EdgeId>     
class PathSet{
public:

    typedef vector<EdgeId> Path ;

    EdgeId start; //These "start" and "end" maybe redundant at this step but will be useful when no-paths are found.
    EdgeId end;
    double length;
    set<Path> paths;

    PathSet(EdgeId start, EdgeId end, double length, set<Path> paths):
        start(start), end(end), length(length), paths(paths)
    {}
    bool operator<(const PathSet& rhs) const {
        const PathSet &lhs = *this;
        if(lhs.start == rhs.start)
        {
            if(lhs.end == rhs.end)
            {
                if(lhs.length == rhs.length)
                {
                    return PathLess(lhs.paths, rhs.paths);
                }
                else return lhs.length< rhs.length;
            }
            else return lhs.end < rhs.end;
        }
        else return lhs.start < rhs.start;
    }
private:

    bool PathLess(set<Path> firstSet, set<Path> secondSet)
    {
        //TODO use the lexicography default comparision.
        if(firstSet.size() < secondSet.size())
            return true;
        vector<EdgeId> v1;
        for(auto iter = firstSet.begin() ; iter != firstSet.end() ; ++iter)
        {
            v1.push_back(*iter);
        }
        vector<EdgeId> v2;
        for(auto iter = firstSet.begin() ; iter != firstSet.end() ; ++iter)
        {
            v2.push_back(*iter);
        }
        if(v1.size() < v2.size())
            return true;
        if(v1.size() > v2.size())
            return false;
        for(size_t i = 0 ; i < v1.size()  ; ++i)
        {
            if(v1[i] < v2[i])
                return true;
        }
        return false;
    }
    bool operator==(const PathSet& rhs) const {
		const PathSet &lhs = *this;
		return !(lhs < rhs || rhs < lhs);
	}
    bool operator!=(const PathSet& rhs) const {
		return !(*this == rhs);
	}

};

template<typename EdgeId>
ostream& operator<<(ostream& os, const PathSet<EdgeId>& pathSet) {
    stringstream pathsString;
    for(int i = 0 ; i < pathSet.paths.size() ; ++i)
    {
        pathsString << endl;
        for(int j = 0 ; j < pathSet.paths[i].size() ; j++)
        {
            pathsString << pathSet.paths[i][j] << " " ;
        }
        pathsString<<endl;
    }
    return os << "Start = " << pathSet.start << ", End = " << pathSet.end<< pathsString.str() ;
}

template<typename EdgeId>
const PathSet<EdgeId> MinPathSet(EdgeId id) {

    set<vector<EdgeId> > paths;
	return PathSet<EdgeId>(id, (EdgeId) 0/*numeric_limits<EdgeId>::min()*/,
			numeric_limits<double>::min(), paths);
}

template<typename EdgeId>
const PathSet<EdgeId> MaxPathSet(EdgeId id) {

    set<vector<EdgeId> > paths;
	return PathSet<EdgeId>(id, (EdgeId) -1/*numeric_limits<EdgeId>::max()*/,
			numeric_limits<double>::max(), paths);
}
template<typename EdgeId>
const PathSet<EdgeId> MinPathSet(EdgeId e1, EdgeId e2) {
	PathSet<EdgeId> pathset = MinPathSet(e1);
	pathset.end = e2;
	return pathset;
}
template<typename EdgeId>
const PathSet<EdgeId> MaxPathSet(EdgeId e1, EdgeId e2) {
	PathSet<EdgeId> pathset = MaxPathSet(e1);
	pathset.end = e2;
	return pathset;
}

//Modified from PairInfoIndexData

template<typename EdgeId>
class PathSetIndexData {
public:
	typedef set<PathSet<EdgeId>> Data;
	typedef typename Data::iterator data_iterator;
	typedef typename Data::const_iterator data_const_iterator;
	typedef vector<PathSet<EdgeId>> PathSets;

	typedef std::pair<data_const_iterator, data_const_iterator> iterator_range;

	data_iterator begin() const {
		return data_.begin();
	}

	data_iterator end() const {
		return data_.end();
	}

	size_t size() const {
		return data_.size();
	}

	void AddPathSet(const PathSet<EdgeId>& pathSet) {
		data_.insert(pathSet);
	}

	void DeletePathSet(PathSet<EdgeId>& pathSet) {
        data_.erase(pathSet);
	}

	PathSets GetPathSets(EdgeId e) const {
		return PathSets(LowerBound(e), UpperBound(e));
	}

	PathSets GetPathSets(EdgeId e1, EdgeId e2) const {
		return PathSets(LowerBound(e1, e2), UpperBound(e1, e2));
	}


	void clear() {
		data_.clear();
	}

	data_iterator LowerBound(EdgeId e) const {
		return data_.lower_bound(MinPathSet(e));
	}

	data_iterator UpperBound(EdgeId e) const {
		return data_.upper_bound(MaxPathSet(e));
	}

	data_iterator LowerBound(EdgeId e1, EdgeId e2) const {
		return data_.lower_bound(MinPathSet(e1, e2)); 
	}

	data_iterator UpperBound(EdgeId e1, EdgeId e2) const {
		return data_.upper_bound(MaxPathSet(e1, e2));
	}

private:
	Data data_;
};

}
