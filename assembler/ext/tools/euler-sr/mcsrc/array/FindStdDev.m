function [mu, sigmaLow, sigmaHigh] = FindStdDev(cumdist)
% FindStdDev computes the standard deviation based on the cumulative
%   distribution of a more or less normal sample.  
  nbins = size(cumdist,2);
  lowIndex = 0;
  meanIndex = 0;
  highIndex = 0;
  for i=1:nbins-1
%    fprintf('cdi: %f  cdi+1: %f\n', cumdist(i), cumdist(i+1));
    if (cumdist(i) <= 0.158 && cumdist(i+1) > 0.158)
      lowIndex = i;
    end
    if (cumdist(i) <= 0.50 && cumdist(i+1) > 0.50)
      meanIndex = i;
    end
    if (cumdist(i) <= 0.80 && cumdist(i+1) > 0.80)
      highIndex = i;
    end
  end
  sigmaLow = 0;
  sigmaHigh = 0;
  mu = 0;
  if (lowIndex > 0 && meanIndex > 0)
    sigmaLow = meanIndex - lowIndex + 1;
  end
  
  if (highIndex > 0 && meanIndex > 0)
    sigmaHigh = highIndex -  meanIndex  + 1;
  end
  fprintf('%d %d %d\n', sigmaLow, meanIndex, sigmaHigh);
  mu = meanIndex;
  
  
    
    
    
    
