
template<class Graph>
class PathSetStats {

	typedef typename Graph::EdgeId EdgeId;
	const Graph& g_;
    const PathSetIndexData<EdgeId> &rawData_;
    const PathSetIndexData<EdgeId> &filteredData_;
private:
void PathsNumbersStats()
{
    vector<size_t> rawVector;
    vector<size_t> filteredVector;

    for(auto iter = rawData_.begin() ; iter != rawData_.end() ; ++iter)
    {
        PathSet<EdgeId> pathset = *iter;
        rawVector.push_back( pathset.paths.size());
    }
    for(auto iter = filteredData_.begin() ; iter != filteredData_.end() ; ++iter)
    {
        PathSet<EdgeId> pathset = *iter;
        filteredVector.push_back(pathset.paths.size());
    }
    //any function to write these vector to files for plotting hist?
}
public:

	PathSetStats(const Graph& g,const PathSetIndexData<EdgeId> &rawData, const PathSetIndexData<EdgeId> &filteredData):
        g_(g), rawData_(rawData), filteredData_(filteredData)
    {
    }

	void Count() {
        PathsNumbersStats();
    }


};
