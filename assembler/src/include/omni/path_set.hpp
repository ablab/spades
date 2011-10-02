#pragma once
#include <cmath>
#include <map>
#include <xmath.h>

namespace omnigraph{

/**
 * PathSet: the basic data structure for the path-set graph
 **/
template<typename EdgeId>     
class PathSet{
public:

    typedef vector<EdgeId> Path ;

    EdgeId start;
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
                if(lhs.distance == rhs.distance)
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


}
