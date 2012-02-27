function dist = DistanceMatrix(mat) 
  nrows = size(mat,1);
  ncols = size(mat,2);
  dist = zeros(nrows*(nrows-1)/2,1);
  isnanval = isnan(mat);
  numnan   = sum(isnan(mat'));
  tic
  index = 1;
  for i = 1:nrows
  if (mod(i,100) == 0)
      toc
      fprintf('i: %d\n', i);
      tic
    end
    for j = i+1:nrows
      if (numFinite > 0)
	meanDist = measuredDist / numFinite;
	dist(index) = measuredDist + meanDist * numNan;
      else
	dist(index) = NaN;
      end
    end
    index = index + 1;
  end
      
      
