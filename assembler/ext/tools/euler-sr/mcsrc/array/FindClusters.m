function FindClusters(tree, index, distrib, distmean, distvar,  nexp, nboot)
  nvertex = size(tree,1);
  if (tree(index,1) > nvertex && tree(index,2) > nvertex)
    % found two siblings in the tree.
    clone1 = tree(index,1);
    clone2 = tree(index,2);
    sampsize = size(distrib,3);
    for e = 1:nexp
      % grab the two experiments
      mu1 = distmean(clone1, e);
      mu2 = distmean(clone2, e);
      sigma12 = distvar(clone1, e);
      if (~isnan(sigma12))
	sigma1 = sqrt(sigma1);
      end
      sigma22 = distvar(clone2, e);
      if (~isnan(sigma22))
	sigma2 = sqrt(sigma2);
      end
      if (~ isnan(mu1) && ~ isnan(mu2) )
	exp1(1:sampsize) = distrib(clone1, e, 1:sampsize);
	exp2(1:sampsize) = distrib(clone2, e, 1:sampsize);
	T = (mu1 - mu2) / ( (sigma12 + sigma22) / 2)  )
      end
    end
  end
      
      
