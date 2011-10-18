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
    int id;

    PathSet(EdgeId start, EdgeId end, double length, set<Path> paths):
        start(start), end(end), length(length), paths(paths), id(-1)
    {}


    bool operator<(const PathSet& rhs) const {
        const PathSet &lhs = *this;
        
        
        
        if(lhs.start == rhs.start)
        {
            if(lhs.length == rhs.length)
            {
                if(lhs.end== rhs.end)
                {
                    return PathLess(lhs.paths, rhs.paths);
                }
                else return lhs.end< rhs.end;
            }
            else return lhs.length < rhs.length;
        }
        else return lhs.start < rhs.start;




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

    bool IsAbsolutePrefixOf(const PathSet<EdgeId> rhs)
    {
        const PathSet<EdgeId> &lhs = *this;
        if(rhs.start != lhs.start)
            return false;
        if(lhs.length >= rhs.length)
            return false;
        for(auto iter = lhs.paths.begin() ; iter != lhs.paths.end() ; ++iter)
        {
            Path singlePath( iter->begin() , iter->end());
            singlePath.push_back(end);
            if (!IsPrefixOfPaths(singlePath, rhs.paths))
                return false ;
        }
        return true;
    }
    bool IsAbsoluteSuffixOf(const PathSet<EdgeId> rhs)
    {
        return true;
    }
    bool CanBeExtendedBy(const PathSet<EdgeId> rhs)
    {
        return true;
    }
    


private:
    bool IsPrefixOfPaths(Path & singlePath,const set<Path> &paths)
    {
        for(auto iter = paths.begin() ; iter != paths.end(); ++iter)
        {
            if( singlePath.size() > iter->size())
                continue;
            bool diff = true;
            for(size_t i = 0  ; i < singlePath.size() ; i++)
            {
                if(singlePath[i] != (*iter)[i])
                {
                    diff = false;
                    break;
                }
            }
            if(diff == true)
                return true;
        }
        return false;
    }


    bool PathLess(const set<Path> firstSet,const set<Path> secondSet) const
    {
        //TODO use the lexicography default comparision.
        if(firstSet.size() < secondSet.size())
            return true;
        if(firstSet.size() > secondSet.size())
            return false;

        vector<Path> v1;
        for(auto iter = firstSet.begin() ; iter != firstSet.end() ; ++iter)
        {
            v1.push_back(*iter);
        }
        vector<Path> v2;
        for(auto iter = firstSet.begin() ; iter != firstSet.end() ; ++iter)
        {
            v2.push_back(*iter);
        }

        for(size_t i = 0 ; i < v1.size()  ; ++i)
        {
            if(v1[i] < v2[i]) //lexicalgraphic comparison
                return true;
        }
        return false;
    }
public:

    bool operator==(const PathSet& rhs) const {
		const PathSet &lhs = *this;
		return !(lhs < rhs || rhs < lhs);
	}
    bool operator!=(const PathSet& rhs) const {
		return !(*this == rhs);
	}
    void SetId(int id_to_set) {
    	id = id_to_set;
    }

};

template<typename EdgeId>
ostream& operator<<(ostream& os, const PathSet<EdgeId>& pathSet) {
    stringstream pathsString;
    size_t linecounter = 1;

    for(auto iter = pathSet.paths.begin() ; iter != pathSet.paths.end() ; ++iter)
    {
        pathsString << "Path " << linecounter <<":"<< pathSet.length<< " "<<  pathSet.start <<"--" ;
        for(size_t i = 0 ; i < (*iter).size() ; ++i)
        {
            pathsString << (*iter)[i] << " -- " ;
        }
        pathsString<<  pathSet.end;
        pathsString<<endl;
    }

    return os << "id: "<< pathSet.id <<" Start = " << pathSet.start <<" ....... "<<"End = " << pathSet.end<< endl<< pathsString.str() ;
}

template<typename EdgeId>
const PathSet<EdgeId> MinPathSet(EdgeId eid) {

    set<vector<EdgeId> > paths;
	return PathSet<EdgeId>(eid, (EdgeId) 0/*numeric_limits<EdgeId>::min()*/,
			numeric_limits<double>::min(), paths);
}

template<typename EdgeId>
const PathSet<EdgeId> MaxPathSet(EdgeId eid) {

    set<vector<EdgeId> > paths;
	return PathSet<EdgeId>(eid, (EdgeId) -1/*numeric_limits<EdgeId>::max()*/,
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
    PathSetIndexData(){
    	maxId = 0;
    	data_.clear();
    }

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

	void AddPathSet(PathSet<EdgeId>& pathSet) {
		pathSet.SetId(maxId);
		maxId++;
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
	int maxId;
};
template<typename EdgeId>
class PathSetIndex
{
private:
    PathSetIndexData<EdgeId> data_;
    typedef vector<EdgeId> Path;

public:
    PathSetIndex(PathSetIndexData<EdgeId> data):data_(data){}
    //TODO BAD CODE

    
    void Process(PathSetIndexData<EdgeId>& filteredPathSetData)
    {
        PathSetIndexData<EdgeId> removedPrefixData;
        RemovePrefixes(removedPrefixData);
        RemoveInvalidPaths(removedPrefixData, filteredPathSetData);
    }


    /*
     * Find a vector of possible extension pathsets.
     * The offSet will take care of the case where there is no pairedinfo in a very short edge. I ignore the offset by now
     * and will consider later when we have problem with the pairinfo. 
     * ****** Seems that we can not ignore the pathset in reality. There are something that has no pairinfo
     */
    void FindExtension(PathSetIndexData<EdgeId> &indexData, PathSet<EdgeId> &currentPathSet, vector<PathSet<EdgeId>> & extendablePathSets, size_t offSet = 0)
    {

        if(currentPathSet.paths.size() == 1 && (*(currentPathSet.paths.begin())).size() ==0 )
        {
            if(offSet > 0)
                return;

            EdgeId headEdge = currentPathSet.end;
            vector<PathSet<EdgeId>>  pathSets =  indexData.GetPathSets(headEdge);
            extendablePathSets.reserve(pathSets.size());
            copy(pathSets.begin(), pathSets.end(), back_inserter(extendablePathSets));

        }
        else
        {
            for(auto iter = currentPathSet.paths.begin() ; iter != currentPathSet.paths.end() ; ++iter)
            {
                Path currentPath  = *iter;
                if(currentPath.size() == 0)
                {
                    INFO("NEED TO BE INSPECTED");
                }
                else{
                    currentPath.push_back( currentPathSet.end);
                    if(currentPath.size() < offSet + 1)
                        continue;
                    Path comparedPath(currentPath.begin() +1 + offSet, currentPath.end());
                    vector<PathSet<EdgeId>> candidatePathSets = indexData.GetPathSets(currentPath[offSet]);
                    for(auto it = candidatePathSets.begin() ; it != candidatePathSets.end() ; ++it)
                    {
                        set<Path> comparedPathsCandidate;
                        for(auto piter  = it->paths.begin() ; piter != it->paths.end() ; ++piter)
                        {
                            Path extendedPath(piter->begin(), piter->end());
                            extendedPath.push_back(it->end);
                            comparedPathsCandidate.insert(extendedPath);
                        }
                        if(IsPrefixOfPaths(comparedPath, comparedPathsCandidate))
                            if(find(extendablePathSets.begin(), extendablePathSets.end(), *it) == extendablePathSets.end()){
                                extendablePathSets.push_back(*it);
                            }
                    }
                }
            }
        
        }
    }

    void RemovePrefixes(PathSetIndexData<EdgeId> &filteredPathSetDat)
    {
        for(auto iter = data_.begin() ; iter != data_.end() ; )
        {
            int distance =0;
            PathSet<EdgeId> currentPathset = *iter;
            auto forward_iter = iter ;  
            bool isPrefix = false ;
            while(true)
            {
                if(iter == data_.end())
                    break;

                advance(iter,1);
                distance++;
            
                if((iter == data_.end()) || ( iter->start != currentPathset.start))
                    break;
                if(currentPathset.IsAbsolutePrefixOf( *iter))
                {
                    isPrefix = true;
                    break;
                }
            }
           if(!isPrefix)
           {
                filteredPathSetDat.AddPathSet(currentPathset);
           }

            advance(iter, -1*distance + 1 );

        }
    }
private:

    void RemoveInvalidPaths(PathSetIndexData<EdgeId> &rawPathSetDat, PathSetIndexData<EdgeId> &filtered)
    {
        //we use ad hoc paired de Bruijn graph for mate-pair info just in the paper
        //Here we use the path-set platform too
        for(auto iter = rawPathSetDat.begin() ; iter != rawPathSetDat.end() ;++iter)
        {
            PathSet<EdgeId> rawPathSet = *iter;
            PathSet<EdgeId> newPathSet = *iter;
            if(rawPathSet.paths.size() == 1)
            {
                filtered.AddPathSet(newPathSet);
                continue;
            }
            vector<PathSet<EdgeId>> topLevelNodes ;
            topLevelNodes.push_back(rawPathSet);
            for(auto pathIter = rawPathSet.paths.begin() ; pathIter != rawPathSet.paths.end() ; ++pathIter)
            {
                Path currentPath = *pathIter;
                deque<EdgeId> checkPath(pathIter->begin(), pathIter->end());
                checkPath.push_front(rawPathSet.start);
                checkPath.push_back(rawPathSet.end);
                if(checkPath.size() == 2)
                    continue;
                if(!CheckForwardConsistent(checkPath, topLevelNodes, rawPathSetDat))
                {
                    newPathSet.paths.erase(currentPath);
                }
            }
            if(newPathSet.paths.size() ==0)
            {
                INFO("ALL PATHS IS REMOVED ---- ");
               // filtered.AddPathSet(rawPathSet);
                INFO(rawPathSet);
            }
            else
            {
                INFO("PATHSET IS VALID ");
                INFO("BEFORE:");
                INFO(rawPathSet);
                INFO("AFTER:");
                filtered.AddPathSet(newPathSet);
                INFO(newPathSet);
            }
        }
    }
    bool CheckForwardConsistent(deque<EdgeId> checkPath, vector<PathSet<EdgeId>> currentPathsets, PathSetIndexData<EdgeId> &pathsetData)
    {
        //base case
        if(checkPath.size() == 1)  
        {
            for(size_t i = 0 ; i < currentPathsets.size() ; ++i )
            {
                if(currentPathsets[i].start == checkPath[0])
                    return true;
            }
            return false;
        }
        else
        {
            EdgeId headNode = checkPath[0];
            vector<PathSet<EdgeId>> nextLevelPathSets ;
            
            size_t offSet = 0;
            vector<PathSet<EdgeId>>  allpossiblePathSets =  pathsetData.GetPathSets(headNode);
            if((allpossiblePathSets.size() == 0) && checkPath.size() > 2)
                offSet = 1;

            for(size_t i =  0 ; i < currentPathsets.size() ; ++i)
            {
                if(currentPathsets[i].start == headNode)
                {
                    //if this is actually a very short edge and there is no pairInfo, we can skip to the next, but not more
   //             EdgeId headEdge = currentPathSet.end;
   

                   vector<PathSet<EdgeId>> extendablePathSets; 
                    FindExtension(pathsetData, currentPathsets[i], extendablePathSets,offSet);
                    for(size_t t = 0 ; t< extendablePathSets.size() ; t++)
                    {
                        if(find(nextLevelPathSets.begin(), nextLevelPathSets.end(), extendablePathSets[t]) == nextLevelPathSets.end())
                            nextLevelPathSets.push_back(extendablePathSets[t]);    
                    }
                }
            }
            checkPath.pop_front();
            if(checkPath.size() > 1 && offSet == 1)
                checkPath.pop_front();
            if(nextLevelPathSets.size() == 0)
                return false;
            else
            {
                if(checkPath.size() == 1)
                    return true;
                return CheckForwardConsistent(checkPath, nextLevelPathSets , pathsetData);
            
            }

        }

    }

   //TODO Move out duplicated functions --- lazy coder 
    
    bool IsPrefixOfPaths(Path & singlePath,const set<Path> &paths)
    {
        for(auto iter = paths.begin() ; iter != paths.end(); ++iter)
        {
            if( singlePath.size() > iter->size())
                continue;
            bool diff = true;
            for(size_t i = 0  ; i < singlePath.size() ; i++)
            {
                if(singlePath[i] != (*iter)[i])
                {
                    diff = false;
                    break;
                }
            }
            if(diff == true)
                return true;
        }
        return false;
    }
    
};

}
