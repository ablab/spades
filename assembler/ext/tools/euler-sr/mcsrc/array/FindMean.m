function  meanvals = FindMean(array, nexp)
  isn = 1 - isnan(array);
  nclones = size(array,1);
  expsize = size(array,2) / nexp;
  for i = 1:nclones
    for j = 1:nexp
      sumV = 0;
      numNonNaN = sum(isn(i,(j-1)*expsize + 1 : j*expsize));
      for k = 1:expsize
	if (isn(i, (j-1)*expsize + k))
	  sumV = sumV + array(i, (j-1)*expsize +k);
	end
      end
      if (numNonNaN == 0)
	meanvals(i,j) = NaN;
      else
	meanvals(i,j) = sumV / numNonNaN;
      end
    end
  end
  
