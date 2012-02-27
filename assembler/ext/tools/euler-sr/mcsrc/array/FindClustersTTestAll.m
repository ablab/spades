function FindClustersTTest(tree, index, expdef, expvals,names)
  nvertex = size(tree,1);
  nexp    = size(expdef,1);
  expSize(1:nexp) = expdef(:,2) - expdef(:,1) + 1;
  if (tree(index,1) < nvertex && tree(index,2) < nvertex)
    % found two siblings in the tree. 
    clone1 = tree(index,1);
    clone2 = tree(index,2);
%    sampsize = size(distrib,3);
%    nexp = size(expdef,1);
    fprintf('%d checking leaves: %d %d clones %s vs %s \n\n', index, ...
	    tree(index,1), tree(index,2), names{clone1}, names{clone2});
    nEqual = 0;
    clear pvalues;
    for expindex = 1:nexp
      % grab the two experiments
      exp1 = SelectExperiment(expvals, clone1, expindex, expdef);
      exp2 = SelectExperiment(expvals, clone2, expindex, expdef);
      sizeExp1 = size(exp1,2);
      sizeExp2 = size(exp2,2);
      if ((sizeExp1 == expSize(expindex)) && (sizeExp2 == ...
					      expSize(expindex)))
%	exp1
%	exp2
	mu1 = mean(exp1);
	mu2 = mean(exp2);
	[rejectEqualMean, pValue] = ttest2(exp1, exp2);
	fprintf('ttest result: %d %f mu1: %f mu2: %f\n', rejectEqualMean, pValue, ...
		mu1, mu2);
	nEqual = nEqual + rejectEqualMean;
	pvalues(size(pvalues,1)+1) = pValue;
      end
    end
    % some metric required here for deciding if the means are equal or not
  else
    if (tree(index,1) > nvertex)
      fprintf('checking children: at %d to %d\n', index, tree(index,1)- ...
	      nvertex-1);
      if (tree(index,1) - nvertex - 1 > 0)
	FindClustersTTest(tree, tree(index,1) - nvertex - 1, expdef, ...
			  expvals, names)
      end
    end
    if (tree(index,2) > nvertex)
      fprintf('checking children: at %d to %d\n', index, tree(index,2) - ...
	      nvertex - 1);
      if (tree(index,2) - nvertex - 1 > 0)
	FindClustersTTest(tree, tree(index,2) - nvertex - 1, expdef, ...
			  expvals, names)
      end
    end
  end

      
      
