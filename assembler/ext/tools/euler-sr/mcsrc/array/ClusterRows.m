function clustering = ClusterRows(data) 
% CLUSTERROWS  Cluster rows (very slowly) based on aggregate distance. 
%
%  Input:
%     clustering = CLUSTERROWS(data) Compute the clustering of
%     distances.  data is a nxp matrix of values to be clustered.
%     
  nrows = size(data,1);
  ncols = size(data,2);
  nClusters = nrows;

  % Each point is in its own cluster, and is centered at itself.  Like me.
  clusterSizes   = ones(1, nrows); % the size of every cluster
  clusterCenters = data;          % the center of each cluster
  clusterIndices = 1:nrows;       % the index of each cluster, to be
                                  % stored by the parent cluster
  clusters       = repmat([0 0 1 1],nrows, 1); % the index of the cluster of
                                          % each (0) for leaves, and the
                                          % distance, all equal to 1 for now.
  
  % merge clusters until there is just one
  while (nClusters > 1) 
    % For each cluster, find its closest neighbor.
    minDist      = repmat(Inf, nClusters, 1);
    maxDist      = repmat(0,   nClusters, 1);
    minDistIndex = repmat(-1, nClusters, 1);
    for i = 1:nClusters-1
      for j = i+1:nClusters
	dist = VectorDistanceWNaN(clusterCenters(i,:), clusterCenters(j,:));
	if (dist < minDist(i)) 
	  minDist(i) = dist;
	  minDistIndex(i) = j;
	end
	% if the min dist index for j has not been assigned
	% assume it is closest to this i
	if (dist < minDist(j) )
	  minDistIndex(j) = i;
	  minDist(j) = dist;
	end
	if (dist > maxDist(i))
	  maxDist(i) = dist;
	end
	if (dist > maxDist(j))
	  maxDist(j) = dist;
	end
      end
    end


    % now when the neighbors are distinct (a->b and no other c->b), 
    % place them on a list to join together.
    numJoined = 0
    clear toJoin;
    for i = 1:nClusters-1
      minIndexI = minDistIndex(i);
      % the two clusters should be merged
      if (minIndexI >= 1 && minDistIndex(minIndexI) == i && minIndexI > i)
	numJoined = numJoined + 1;
	toJoin(numJoined,1:2) = [i, minDistIndex(i)];
      end
    end
    if (numJoined == 0)
      break;
    end
    fprintf('tojoin:\n');
    toJoin
    % now join together clusters
    i = 1;
    joinIndex = 1;
    clusterIndex = size(clusters,1);
    newClusterIndex = 0;
    for i = 1:nClusters
      if (joinIndex == numJoined || i ~= toJoin(joinIndex,1))
	% don't join this cluster, just append it to the list of new
        % clusters.
	if (size(find(toJoin(:,2)==i),1) == 0)
	  newClusterIndex                      = newClusterIndex + 1;
	  newClusterSizes(newClusterIndex)     = clusterSizes(i);
	  newClusterCenters(newClusterIndex,1:ncols) = clusterCenters(i,:);
	  newClusterIndices(newClusterIndex)   = clusterIndices(i);
	else
	  fprintf('skipping %i\n', i);
	end
      else
	% create a new cluster 
	fprintf('some debugging information i: %d ji: %d nj: %d\n', i, ...
		joinIndex, numJoined);
	j = toJoin(joinIndex,2);
	newClusterIndex = newClusterIndex + 1;
	newClusterSizes(newClusterIndex) = clusterSizes(i) + ...
	    clusterSizes(j);
	newClusterCenters(newClusterIndex,1:ncols) = ...
	    AddWeightedVectors(clusterCenters(i,:), ...
			       clusterSizes(i), ...
			       clusterCenters(j,:),...
			       clusterSizes(j));

%	fprintf('created new cluster center\n');
%	newClusterCenters(newClusterIndex,1:ncols)
%	fprintf('from the following two:\n');
%	clusterCenters(i,:)
%	clusterCenters(j,:)
	clusterIndex = clusterIndex + 1;
	newClusterIndices(newClusterIndex) = clusterIndex;
	% add this cluster to the list of clusters
	clusters(clusterIndex,1:4) = [clusterIndices(i) clusterIndices(j) ...
		    minDist(i) maxDist(i)];
	joinIndex = joinIndex + 1;
      end
    end
    fprintf('done, swapping values\n');
    % swap the old data out.
    oldcs = size(clusterSizes,2);
    oldcc = size(clusterCenters,1);
    oldci = size(clusterIndices,2);
    clear clusterSizes;
    clear custerCenters;
    clear clusterIndices;
    clusterSizes   = newClusterSizes;
    clusterCenters = newClusterCenters;
    clusterIndices = newClusterIndices;
    clear newClusterSizes;
    clear newClusterCenters;
    clear newClusterIndices;
    newcs = size(clusterSizes,2);
    newcc = size(clusterCenters,1);
    newci = size(clusterIndices,2);
    fprintf('old: cs: %d cc: %d ci: %d  new: cs: %d cc: %d ci: %d\n', ...
	    oldcs, oldcc, oldci, newcs, newcc, newci);
    fprintf('new cluster index: %d\n', newClusterIndex);
    nClusters = size(clusterSizes,2)
  end
  totalClust = size(clusters,1);
  size(clusters)
  clustering = clusters(nrows+1:totalClust,1:3);

  
      
    
    
	    
	
