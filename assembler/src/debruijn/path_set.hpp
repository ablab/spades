#pragma once
#include <cmath>
#include <map>
#include <limits>
#include <xmath.h>
#include "graph_pack.hpp"
//#include "path_set_tools.hpp"
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
    double weight;
    PathSet(EdgeId start, EdgeId end, double length, set<Path> paths, double weight):
        start(start), end(end), length(length), paths(paths), id(-1), weight(weight)
    {}

    template<class Comparator>
    class PathSetComparator {
    private:
    	Comparator comparator_;
    public:
    	PathSetComparator(Comparator comparator) : comparator_(comparator) {
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

    	bool operator()(const PathSet& lhs, const PathSet& rhs) {
            if(lhs.start == rhs.start)
            {
                if(lhs.length == rhs.length)
                {
                    if(lhs.end== rhs.end)
                    {
                        return PathLess(lhs.paths, rhs.paths);
                    }
                    else return comparator_(lhs.end, rhs.end);
                }
                else return comparator_(lhs.length, rhs.length);
            }
            else return comparator_(lhs.start, rhs.start);




            if(lhs.start == rhs.start)
            {
                if(lhs.end == rhs.end)
                {
                    if(lhs.length == rhs.length)
                    {
                        return PathLess(lhs.paths, rhs.paths);
                    }
                    else return comparator_(lhs.length, rhs.length);
                }
                else return comparator_(lhs.end, rhs.end);
            }
            else return comparator_(lhs.start, rhs.start);
    	}
    };

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
/*
template<typename EdgeId>
ostream& operator<<(ostream& os, const PathSet<EdgeId>& pathSet) {
    stringstream pathsString;
    size_t linecounter = 1;

    for(auto iter = pathSet.paths.begin() ; iter != pathSet.paths.end() ; ++iter)
    {
        pathsString << "Path " << linecounter <<":"<< pathSet.length<< " "<<  pathSet.start <<"--" ;
        linecounter++;
        for(size_t i = 0 ; i < (*iter).size() ; ++i)
        {
            pathsString << (*iter)[i] << " -- " ;
        }
        pathsString<<  pathSet.end;
        pathsString<<endl;
    }

    return os << "id: "<< pathSet.id <<" weight "<< pathSet.weight <<" Start = " << pathSet.start <<" ....... "<<"End = " << pathSet.end<< endl<< pathsString.str() ;
}
*/
template<typename EdgeId>
const PathSet<EdgeId> MinPathSet(EdgeId eid) {

    set<vector<EdgeId> > paths;
	return PathSet<EdgeId>(eid, (EdgeId) 0/*numeric_limits<EdgeId>::min()*/,
			numeric_limits<double>::min(), paths, 0);
}

template<typename EdgeId>
const PathSet<EdgeId> MaxPathSet(EdgeId eid) {

    set<vector<EdgeId> > paths;
	return PathSet<EdgeId>(eid, (EdgeId) -1/*numeric_limits<EdgeId>::max()*/,
			numeric_limits<double>::max(), paths, 0);
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

	int AddPathSet(PathSet<EdgeId>& pathSet) {
		pathSet.SetId(maxId);
		int tmp = maxId;
		maxId++;
		data_.insert(pathSet);
		return tmp;
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


template<typename graph_pack>
class PathSetIndex
{
	typedef typename graph_pack::graph_t::EdgeId EdgeId;
	typedef typename graph_pack::graph_t::VertexId VertexId;
	typedef typename graph_pack::graph_t Graph;
private:
    PathSetIndexData<EdgeId> data_;
    PathSetIndexData<EdgeId> backwarddata_;
    vector<PathSetIndexData<EdgeId> > additional_libs;
    typedef vector<EdgeId> Path;
    graph_pack& gp;
    
    void GenerateBackwardData()
    {
        PathSetIndexData<EdgeId> rawBackwardData;
        for(auto iter = data_.begin() ; iter != data_.end() ; ++iter)
        {
            PathSet<EdgeId> forwardPathSet = *iter;
            set<Path> paths ;
            for(auto itPath = forwardPathSet.paths.begin() ; itPath != forwardPathSet.paths.end() ; ++itPath)
            {
                Path apath = *itPath;
                reverse(apath.begin(), apath.end());
                paths.insert(apath);
            }

            PathSet<EdgeId> pathset(forwardPathSet.end , forwardPathSet.start, forwardPathSet.length , paths, forwardPathSet.weight);
            rawBackwardData.AddPathSet(pathset);

        }
        for(auto iter = rawBackwardData.begin() ; iter != rawBackwardData.end() ; )
        {
            int distance =0;
            PathSet<EdgeId> currentPathset = *iter;
//            auto forward_iter = iter ;
            bool isPrefix = false ;
            while(true)
            {
                if(iter == rawBackwardData.end())
                    break;

                advance(iter,1);
                distance++;

                if((iter == rawBackwardData.end()) || ( iter->start != currentPathset.start))
                    break;
                if(currentPathset.IsAbsolutePrefixOf( *iter))
                {
                    isPrefix = true;
                    //                    iter->SetWeight(iter->weight + currentPathset.weight );
                    break;
                }
            }
            if(!isPrefix)
            {
                backwarddata_.AddPathSet(currentPathset);
            }

            advance(iter, -1*distance + 1 );


        }
    }

public:
    PathSetIndex(PathSetIndexData<EdgeId> data, graph_pack& gp):data_(data), gp(gp){
        GenerateBackwardData();
        additional_libs.clear();
    }

    void AddAdditionalLib(PathSetIndexData<EdgeId>& lib){
    	additional_libs.push_back(lib);
    }

    void Process(PathSetIndexData<EdgeId>& filteredPathSetData)
    {
        PathSetIndexData<EdgeId> removedPrefixData;
        RemovePrefixes(removedPrefixData);
        RemoveInvalidPaths(removedPrefixData, filteredPathSetData);
        PathSetIndexData<EdgeId> finnestData;
        SplitPathSet(filteredPathSetData,finnestData);
    }

    void SplitPathSet(PathSetIndexData<EdgeId>&inputPathSetData, PathSetIndexData<EdgeId> &splitPathSetData )
    {
        for(auto iter = inputPathSetData.begin() ; iter != inputPathSetData.end() ;++iter)
        {
            PathSet<EdgeId> currentPathSet = *iter;
            if(currentPathSet.paths.size() == 1)
            {
                splitPathSetData.AddPathSet(currentPathSet);
            }
            else
            {
                vector<set<Path>> resultsPartition; 
                bool canbeSplitted = TrySplitSingletons(currentPathSet, inputPathSetData,  resultsPartition);
                for(size_t i =0 ; i < resultsPartition.size() ; ++i)
                {
                    if(!canbeSplitted)
                    {
                        splitPathSetData.AddPathSet(currentPathSet);
                        break;
                    }
                    PathSet<EdgeId> pathset(currentPathSet.start , currentPathSet.end, currentPathSet.length , resultsPartition[i], currentPathSet.weight);
                    splitPathSetData.AddPathSet(pathset);
                }
            }
        }
    }

    /*
     * We test if all of these paths are mandatory path. A path is a mandatory path if its removal corrupts the covering walk in the graph.*/
    bool TrySplitSingletons(PathSet<EdgeId> &pathset, PathSetIndexData<EdgeId> &pathSetIndex, vector<set<Path>> &resultsPartition)
    {
        set<Path> paths = pathset.paths;
        assert(paths.size() > 0);
        bool canbeSplit = false ;
        set<EdgeId> allLocaEdges;
        //Conjecture: If there are two paths and they can be verified by back and forward, then both of them are valid ..I can't find a way not using the graph.
        if(paths.size() ==2)
        {
            Path path1 = *paths.begin();
            auto iter = paths.begin();
            ++iter;
            Path path2 = *iter;
            set<Path> set1 ; 
            set<Path> set2;
            set1.insert(path1);
            set2.insert(path2);

            resultsPartition.push_back(set1);
            resultsPartition.push_back(set2);
            return true;
        }
        /*
        for(auto iter = paths.begin() ; iter != paths.end() ; ++iter)
        {
            Path path = *iter;
            for(size_t i = 0 ; i < path.size() ; ++i)
            {
                allLocaEdges.insert(path[i]);
            }
        }
        */

        return canbeSplit;
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
//            auto forward_iter = iter ;
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
//                    iter->SetWeight(iter->weight + currentPathset.weight );
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
            
            vector<PathSet<EdgeId>> topLevelBackwardNodes = backwarddata_.GetPathSets(rawPathSet.end);

            for(auto pathIter = rawPathSet.paths.begin() ; pathIter != rawPathSet.paths.end() ; ++pathIter)
            {
                Path currentPath = *pathIter;
                deque<EdgeId> checkPath(pathIter->begin(), pathIter->end());
                checkPath.push_front(rawPathSet.start);
                checkPath.push_back(rawPathSet.end);
                if(checkPath.size() == 2)
                    continue;
                deque<EdgeId> reverseCheckPath = checkPath;
                reverse(reverseCheckPath.begin(), reverseCheckPath.end());
                
                if(topLevelBackwardNodes.size() == 0)
                {

                    if(!CheckForwardConsistent(checkPath, topLevelNodes, rawPathSetDat))
                    {
                        newPathSet.paths.erase(currentPath);
                    }
                }
                else 
                {
                    if((!CheckForwardConsistent(checkPath, topLevelNodes, rawPathSetDat)) || (!CheckForwardConsistent(reverseCheckPath,topLevelBackwardNodes, backwarddata_ ) ) )
                    {
                        newPathSet.paths.erase(currentPath);
                    }
                }
            }


            if(newPathSet.paths.size() ==0)
            {
                INFO("ALL PATHS IS REMOVED ---- ");
               // filtered.AddPathSet(rawPathSet);
                INFO(str(rawPathSet));
            }
            else
            {
                INFO("PATHSET IS VALID ");
                INFO("BEFORE:");
                INFO(str(rawPathSet));
                INFO("AFTER:");
                filtered.AddPathSet(newPathSet);
                INFO(str(newPathSet));
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
            for(size_t forwardPosition = 1 ; forwardPosition < checkPath.size() ; ++forwardPosition)
            {
                if(pathsetData.GetPathSets(checkPath[forwardPosition]).size() != 0)
                    break;
                offSet++;
            }
            //We just allow 2 skips ?
            if(offSet >= 2)
                return false;


            for(size_t popNum = 0 ; popNum <= offSet ; ++popNum)
            {
                checkPath.pop_front();
            }
            if(checkPath.size() == 0)
                return true;

            for(size_t i =  0 ; i < currentPathsets.size() ; ++i)
            {
                if(currentPathsets[i].start == headNode)
                {
                    //if this is actually a very short edge and there is no pairInfo, we can skip to the next, and next if the condition still holds 
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
    string str(const PathSet<EdgeId> pathSet) {
    	stringstream pathsString;
    	size_t linecounter = 1;
    	for(auto iter = pathSet.paths.begin() ; iter != pathSet.paths.end() ; ++iter)
    	{
    		pathsString << "Path " << linecounter <<":"<< pathSet.length<< " "<<  gp.int_ids.ReturnIntId(pathSet.start) <<"--" ;
    		linecounter++;
    		for(size_t i = 0 ; i < (*iter).size() ; ++i)
    		{
    			pathsString << gp.int_ids.ReturnIntId(((*iter)[i])) << " -- " ;
    		}
    		pathsString<<  gp.int_ids.ReturnIntId(pathSet.end);
    		pathsString<<endl;
    	}
    	stringstream res;
    	res << "id: "<< pathSet.id <<" weight "<< pathSet.weight <<" Start = " << gp.int_ids.ReturnIntId(pathSet.start) <<" ....... "<<"End = " << gp.int_ids.ReturnIntId(pathSet.end)<< endl<< pathsString.str() ;
    	return res.str();
    }


    
};

}
